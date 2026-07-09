To structure the new **Photon Package Structure of Arrays (SoA)** for GPU offloading, you should decompose the 11 specific fields currently stored in Enzo’s photon package objects into individual contiguous arrays. This transition from a linked list to an SoA is critical for **coalesced memory access** and to expose the "massive parallelism" required for many-core architectures 1, 2\.  
Based on the specifications of the Moray ray tracer, the SoA should be structured as follows:

### 1\. Defined Data Fields

According to the source documentation, each photon package tracks 11 key attributes 3\. In an SoA format, these would be represented as separate arrays:

* **Floating-Point Arrays (Single Precision):**  
* PhotonFlux\[\]: The number of photons or energy in the package 3\.  
* PhotonEnergy\[\]: The energy of the package (often the average above an ionization threshold) 3, 4\.  
* TotalColumnDensity\[\]: The integrated density traversed by the ray 3\.  
* **Integer Arrays:**  
* PhotonType\[\]: Identifies the radiation category (e.g., H I ionizing, He II ionizing, or Lyman-Werner) 3, 5\.  
* HealpixPixel\[\]: The specific pixel number on the HEALPix sphere 3, 6\.  
* HealpixLevel\[\]: The current refinement level of the ray split 3, 7\.  
* **Double-Precision Arrays (Precision-Dependent):**  
* *Note: If Enzo is compiled with double precision for grid/particle positions, these five items MUST use double precision to ensure spatial accuracy 3\.*  
* TimeInterval\[\]: The emission time-interval 3\.  
* EmissionTime\[\]: The time the photon was created 3\.  
* CurrentTime\[\]: The current internal time of the package 3\.  
* Radius\[\]: The distance from the source 3, 7\.  
* SourcePosition\[8\]\[\]: A set of three arrays ($x, y, z$) for the originating source coordinates 3, 7\.

### 2\. Conceptual C++ Structure

A robust implementation for OpenMP offloading would look like this:  
struct PhotonPackageSoA {  
    int numPackages;  
    int capacity; // Pre-allocated size to handle ray splitting

    // Contiguous arrays for GPU mapping  
    float  \*Flux;  
    int    \*Type;  
    float  \*Energy;  
    double \*TimeInterval;  
    double \*EmissionTime;  
    double \*CurrentTime;  
    double \*Radius;  
    float  \*ColumnDensity;  
    int    \*PixelNum;  
    int    \*Level;  
    double \*SourceX, \*SourceY, \*SourceZ;  
};

### 3\. Key Implementation Strategies

* **Contiguous Memory Mapping:** You must use \#pragma omp target data map(to: ...) to transfer these arrays as a single block to the GPU. This replaces the pointer-chasing overhead of the original doubly linked list 1\.  
* **Pre-allocation for Splitting:** Because rays split into four child rays when sampling becomes too coarse (determined by the cell face area and the ray's solid angle $\\Omega\_{ray}$), your SoA must be pre-allocated with enough "headroom" or capacity to handle the increased number of packages without a CPU-GPU re-allocation mid-step 7, 9\.  
* **Vectorization:** This structure allows the compiler to map each ray propagation loop directly to a GPU thread (SIMD/SIMT). This is essential for calculating **ray segment lengths ($dr$)**, **optical depths ($\\tau$)**, and **photo-ionization rates ($k\_{ph}$)** efficiently in parallel 10, 11\.  
* **Memory Footprint:** A single photon package in double precision is approximately **88 bytes** 3\. When structuring your SoA, ensure that the total memory used by these arrays does not exceed the available device VRAM, especially in large simulations where the number of ray segments can be very high 12\.

