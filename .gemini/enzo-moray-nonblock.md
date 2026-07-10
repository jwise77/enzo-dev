To implement a clean and robust **non-blocking MPI communication** version of the Enzo ray tracer, you should follow the parallelization strategy outlined in the original ENZO+MORAY design, which was specifically engineered to overlap communication with ray propagation 1, 2\.  
The primary challenge you likely faced with race conditions in your \#ifdef NONBLOCKING\_RT implementation is the **distributed termination problem**—knowing when every processor is truly finished when rays can split, move, or be deleted dynamically 3, 4\.

### 1\. The Core Non-Blocking Strategy

Based on the sources, the most effective way to structure this in EvolvePhotons.C and its associated routines is a **six-step asynchronous loop** 377–380:

* **Step 1: Handshake for Message Counts.** Before sending photon data, each processor must determine how many rays it will send to every other processor 5\. Use MPI\_Isend to send these counts and MPI\_Irecv to prepare for incoming counts. This ensures the receiving processor knows exactly how many data messages to expect 5\.  
* **Step 2: Contiguous Packing.** As rays reach a grid boundary, do not send them immediately. Move them to a PhotonMoveList and pack them into **contiguous arrays** for communication 5, 6\. This pairs perfectly with the **Structure of Arrays (SoA)** format, allowing you to send large blocks of ray data efficiently 5\.  
* **Step 3: Prepare Data Receives.** Use the counts received in Step 1 to post the necessary MPI\_Irecv calls for the actual photon data buffers 7\.  
* **Step 4: Asynchronous Data Transfer.** Call MPI\_Isend for the grouped photon data. Because Enzo replicated the grid hierarchy metadata on every processor, each rank can independently identify the correct destination rank for every ray 8, 9\.  
* **Step 5: Overlap Computation and Communication.** This is the "heart" of the loop. While waiting for messages to arrive:  
* Continue tracing rays that are **locally available** in your grids 9\.  
* Use MPI\_Testsome to check if any incoming ray messages have arrived 9\.  
* If a message is complete, unpack the rays into the destination grid's queue and process them immediately 9\.  
* **Only call MPI\_Waitsome** if you have zero local rays left to trace, forcing the processor to idle until new work arrives from a neighbor 9\.  
* **Step 6: Global Termination Check.** The loop ends only when all processors have exhausted their workload 4\. To avoid race conditions here, ensure that "workload exhausted" includes both having no local rays *and* having no outstanding MPI requests (sends or receives) 4\.

### 2\. Guidance to Avoid Race Conditions

* **Buffer Persistence:** A common source of race conditions in non-blocking MPI is modifying or deleting a send buffer before the MPI hardware has finished transmitting it. You must maintain a pool of communication buffers and only reuse them once MPI\_Test or MPI\_Wait confirms the request is complete 10\.  
* **Message Ordering:** In a non-blocking environment, the "count" message must always be processed before the "data" message. By "aggressively draining the message stack" at the start of each iteration, you ensure that unexpected handshake messages from other ranks are handled promptly 7\.  
* **Ray Splitting Synchronization:** Because rays can split into four child rays upon entering a new grid if sampling is too coarse ($A\_{cell} / \\Omega\_{ray} \> 3$), your communication logic must account for the fact that a processor might receive 10 rays but generate 40 3, 11\.

### 3\. Software Design Strategy

To address your "messy design" concerns, move away from scattered MPI\_Test checks and adopt a **Task-Based Manager** approach:

1. **Communication Manager:** Create a dedicated class or set of routines to handle the MPI state machine. EvolvePhotons.C should simply call CommManager-\>SendRays() and CommManager-\>ReceiveRays(), hiding the complexity of the MPI\_Request arrays and buffer management.  
2. **State Tracking:** Instead of complex checks, maintain a simple state for each rank (e.g., SENDING\_COUNTS, RECEIVING\_DATA, IDLE).  
3. **Local vs. Remote Queues:** Use a single "Active Ray Queue" for each grid. Incoming MPI rays are simply pushed onto this queue, making the tracer logic identical regardless of whether a ray was local or just arrived from another node 9\.

By utilizing this **aggressive non-blocking approach** used by the Moray tracer, you can significantly reduce the idle time where processors wait for grid boundary synchronization, allowing the ray tracer to scale more effectively to large processor counts 1, 2\.  
