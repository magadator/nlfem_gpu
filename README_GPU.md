# nlfem — GPU element assembly (Phase 1)

## What this actually is

This is **not** a full GPU rewrite of nlfem. It's the first, verified slice
of one: the **cable element** stiffness/mass/internal-force computation now
runs on the GPU in a single batched kernel launch, wired into the existing
`assembleGlobKMat()` loop. Everything else — beam elements, shell elements
(SIMPLETRI/LMITC3/MITC3), the Newton-Raphson nonlinear loop, Newmark-Beta
time integration, PETSc's `MatSetValues`/`KSPSolve`, MPI, ParMETIS domain
decomposition, FSI transfer, VTK output — is completely unchanged.

Read this before running anything on real results.

## Why it's scoped this way

`nlfem` is validated research code (aortic dissection / parachute FSI). A
mechanical, unverified translation of 20,000+ lines of nonlinear FEM math
to CUDA is how you get silently wrong stresses — the code still runs, still
produces plots, and the physics is quietly off. So this phase does three
things, in order, for exactly one element type:

1. Reads the existing CPU implementation (`comp_Ke_Cable` / `comp_Me_Cable`
   in `src/fem/linearTriElement.cc`) and transcribes it faithfully,
   including quirks that look like they might be bugs (see the comment
   block at the top of `src/gpu/cable_kernels.cu` — in particular
   `btb[i*6+j] = BL[j]*BL[j]*emod`, which uses `BL[j]` twice rather than
   `BL[i]*BL[j]`. That's preserved as-is. If it's not intentional, it's your
   call to fix, not something to quietly "correct" during a port).
2. Fixes exactly one thing that the CPU version gets away with but a GPU
   kernel can't: `comp_Me_Cable` reads `x0_t_`/`x1_t_` member scratch that
   was set by the *previous* call to `comp_Ke_Cable` for the same element,
   in the same loop iteration. That's fine in a serial loop, not fine for
   independent parallel threads. The kernel recomputes it explicitly.
3. Ships with a standalone validation harness (`validate_cable.cc`) that
   diffs GPU output against a plain-CPU reference of the same formulas, so
   you can check the kernel before trusting it on your actual mesh.

Beam and shell elements need the same treatment before they're safe to
port. The MITC3 shell in particular shares mutable state (`faceID`, `x0_0_`,
`x1_0_`, `x2_0_`) across ~8 member functions (`calculateStrainsMITC3`,
`comp_KNL1_MITC3`, `comp_KNL2_MITC3`, `compPKST`,
`calculateDeformationGradient`, `PKST2Cauchy`, `calculateVonMisesStress`) —
that has to be untangled into a pure per-element function first. That's the
next phase, not done here.

## Hardware reality check — Quadro K1200

GPU is a 2016-era Maxwell part: ~640 CUDA cores, 4 GB GDDR5, compute
capability 5.0, no tensor cores, no double-precision throughput advantage
(FP64 runs at roughly 1/32 the FP32 rate on this class of card — nlfem does
everything in `double`). Practical implications:

- **Cable elements**: cheap per element, so this is mostly a latency/launch-
  overhead exercise. You'll only see a real wall-clock win if your mesh has
  a lot of cable elements (thousands+) or you're calling assembly many
  times per run (which you are, once per Newton iteration).
- **Shell/MITC3** (when ported): this is where GPU assembly could
  meaningfully help, since it's the expensive element type per-DOF. But the
  K1200's FP64 throughput is modest, so don't expect the kind of speedup
  you'd see on a datacenter card — expect a moderate win on assembly time,
  not a 10x.
- **The linear solve** (still PETSc/CPU in this phase) is very likely your
  actual bottleneck for large meshes, not element assembly. Profile before
  assuming assembly is where the time goes — `TIME_ASSEMBLY=1` in
  `fem.cc` gives you a timer for exactly this.

Bottom line: this phase is a correctness-first foundation, not a
performance claim. Measure on your actual problem before drawing
conclusions about whether GPU assembly is worth extending further on this
card.

## Build

Requires the CUDA toolkit (`nvcc`) in addition to the existing MPI/PETSc/
ParMETIS toolchain. `nvidia-smi` should show your K1200s.

1. Edit `config.mk`:
   - Set `CUDA_DIR` to your CUDA install (`/usr/local/cuda` by default).
   - `CUDA_ARCH` is already set to `sm_50` for the K1200. Confirm with:
     ```
     nvidia-smi --query-gpu=compute_cap --format=csv
     ```
   - Set `USEGPU=0` to build the original all-CPU codebase, unchanged
     (uses `src/gpu/gpu_assembly_stub.cc` instead of the real kernel — this
     is a true no-op path, not a "GPU that does nothing").
2. Build as before:
   ```
   make clean && make
   ```

## Validate before trusting it

```
make validate_cable
./validate_cable
```

This builds a small standalone binary — no PETSc/MPI/ParMETIS needed —
that runs the kernel on four synthetic cable elements (straight, angled,
pre-stretched) and diffs the result against a plain-CPU transcription of
`comp_Ke_Cable`/`comp_Me_Cable`. Expect:

```
max abs diff (Ke, Me, Fvec, tension) across 4 elements: <something like 1e-13>
PASS
```

If it doesn't say `PASS`, do not run this on real data — open an issue for
yourself and fix the kernel first.

For a second, stronger check on your actual mesh: run one assembly pass
with `NLFEM_NO_GPU=1` set (forces CPU path) and again without it, and diff
the resulting `Kmat_p_`/`Mmat_p_` matrices (or just compare the printed
residual norm / displacement after one Newton iteration — they should
match to solver tolerance).

## Run

Unchanged from before:
```
mpirun -np <N> ./nlfem example/input.dat
```

GPU cable assembly is used automatically when:
- the build was done with `USEGPU=1` (default),
- a CUDA device is detected at runtime, and
- the mesh actually has cable elements.

Otherwise it falls back to the original CPU path silently and correctly —
there's no failure mode where you're missing GPU accel and don't know it
except a one-line message if a GPU call throws mid-run (see
`fem::gpuPrecomputeCableElements()` in `fem.cc`), in which case that
assembly pass, and only that pass, falls back to CPU.

To force CPU-only at runtime without recompiling (useful for A/B timing or
debugging a suspected GPU-path bug):
```
NLFEM_NO_GPU=1 mpirun -np <N> ./nlfem example/input.dat
```

## What's NOT done (be clear-eyed about this)

- Beam element (`comp_Ke_BEBeam`/`comp_Me_BEBeam`) — not ported.
- Shell elements (`SIMPLETRI`, `LMITC3`, `MITC3`) — not ported. MITC3 in
  particular needs the stateful-pipeline refactor described above first.
- The linear system solve is still PETSc/CPU. Removing PETSc entirely (as
  originally asked) means replacing `KSPSolve` with a GPU sparse solver
  (cuSPARSE/cuSOLVER or a hand-rolled preconditioned CG) — that's a
  separate, substantial piece of work, and on a 4 GB K1200 you'd want to
  benchmark whether it's actually faster than PETSc+GAMG on your CPU core
  count before committing to it.
- MPI/ParMETIS domain decomposition is untouched. It's somewhat orthogonal
  to a single-GPU port — you'd still want it if you run across multiple
  nodes, though the domain-decomposition boundaries would need to align
  with per-GPU work if you ever go multi-GPU.

## File map

```
src/gpu/gpu_assembly.h        - data structures + function declarations, no CUDA syntax (safe to #include anywhere)
src/gpu/cable_kernels.cu      - the actual CUDA kernel + host launcher (built when USEGPU=1)
src/gpu/gpu_assembly_stub.cc  - CPU no-op stand-in (built when USEGPU=0)
src/gpu/validate_cable.cc     - standalone CPU-vs-GPU validation harness
src/fem/fem.h                 - added: useGPUAssembly_, gpuCable* buffers, gpuPrecomputeCableElements() declaration
src/fem/fem.cc                - added: gpuPrecomputeCableElements() impl, dispatch in assembleGlobKMat()
config.mk                     - added: USEGPU, NVCC, CUDA_DIR, CUDA_ARCH
Makefile                      - added: nvcc build rule, validate_cable target, conditional GPU sources/link flags
```
