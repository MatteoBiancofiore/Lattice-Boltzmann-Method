# D2Q9 Lattice Boltzmann Method (LBM)
## 1. Project Overview

This project is an open-source, high-performance fluid dynamics solver based on the Lattice Boltzmann Method (LBM), developed from the ground up in standard C. The primary objective is to simulate incompressible, isothermal fluid flows within complex geometries with high numerical accuracy and computational efficiency.

The solver is currently in its early active development phase. The primary focus of this initial phase is establishing a mathematically sound and structurally clean baseline implementation. Our ultimate goal is to build an extensible, robust foundation optimized for massive parallelization across high-performance computing (HPC) architectures

---

## 2. Technical and Implementation Details

The implementation adopts the widely studied two-dimensional, nine-velocity discrete lattice model (**D2Q9**), which strikes an optimal balance between physical accuracy and computational complexity for 2D Navier-Stokes approximations.

### Memory Layout and Data Structures
To overcome the severe latency penalties of multidimensional pointer arrays, the domain is laid out in a single, contiguous 1D array allocated directly on the heap.
* **Row-Major Linearization:** Memory indexing follows standard row-major ordering.
* **HPC-Ready Structure:** This contiguous footprint significantly reduces memory fragmentation and directly mirrors the flat memory models required for unified memory spaces.
### Collision Operator
Inter-particle collisions are modeled using the standard Single Relaxation Time (SRT) **Bhatnagar-Gross-Krook (BGK)** approximation. The relaxation towards the Maxwell-Boltzmann local equilibrium distributions is governed by a dimensionless relaxation parameter, directly tied to the kinematic shear viscosity.
The collision step remains strictly local to each lattice node, providing perfect algorithmic independence.

### Boundary Conditions
In the near future we aim to implement in the solver two complementary boundary approaches depending on the physical interface:
* **Standard Half-Way Bounce-Back:** Enforces no-slip wall conditions along top and bottom solid boundaries, as well as on internal arbitrary obstacles (e.g., bluff bodies, cylinders). Populations attempting to enter a solid cell are reflected back along their inverse discrete directions .
* **Zou-He Boundary Scheme:** Applied to open boundaries (channel inlet and outlet). Zou-He evaluates the non-equilibrium bounce-back principle along the surface normal. It dynamically reconstructs missing incoming distribution functions based on the imposed parabolic velocity profile.
  This guarantees conservation of mass and momentum across open domain boundaries.

### Export Pipeline
The solver exports macroscopic field variables into tabular CSV format for validation and post-processing via a dedicated Python pipeline (using NumPy and Matplotlib, currently under development).

---

## 3. Future Work

The planned engineering roadmap comprises three distinct phases:

### Phase I: Shared-Memory Parallelism (OpenMP)
* Exploit multi-core CPU architectures using OpenMP.
* Establish baseline compute metrics, including Millions of Lattice Node Updates Per Second (MLUPS).

### Phase II: Massive Acceleration on Many-Core Architectures (NVIDIA CUDA)
* Map 2D lattice blocks to CUDA thread blocks, using shared memory banks to optimize streaming exchanges.
* Overlap host-to-device asynchronous memory transfers using CUDA streams during ongoing compute cycles.

### Phase III: 3D Extension and Advanced Visualization
* Expand the computational engine from D2Q9 to standard 3D configurations (**D3Q27**) to simulate realistic turbulent structures, vortex shedding, and 3D aerodynamic profiles.
* Integrate high-throughput asynchronous binary I/O to handle high-resolution transient datasets without stalling the computational kernels.
