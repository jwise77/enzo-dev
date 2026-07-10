To implement the handshake for message counts in Enzo’s ray tracer, you must follow the asynchronous, non-blocking communication protocol detailed in the Moray documentation 1, 2\. This handshake is the first critical phase of the parallelization strategy, ensuring that receiving processors know exactly how many data messages to expect from their neighbors before actual photon data is transmitted 3\.

### 1\. Count and Group Outgoing Rays

Before any MPI calls are made, each processor must scan its PhotonMoveList—the list of rays that have reached grid boundaries or need to be moved to child grids—and **group these rays by their destination processor ID** 2, 3\.

### 2\. Determine the Number of Messages ($N\_{mesg}$)

Because MPI messages have size limits, you must determine how many individual data buffers are required for each destination rank.

* **Packet Size:** The documentation suggests a maximum of $10^5$ rays per MPI message ($N\_{max}$) 3\.  
* **Calculation:** For every destination processor, calculate $N\_{mesg}$ by dividing the total number of rays destined for that rank by $N\_{max}$ 3\.

### 3\. Asynchronous Send and Receive of Counts

This is the "handshake" phase that prevents the code from blocking and avoids the race conditions you encountered:

* **The Sender:** Every rank uses **MPI\_Isend** to send its calculated $N\_{mesg}$ to each corresponding destination rank 3\. This send must be non-blocking so the processor can immediately begin other tasks, such as packing the photon data arrays 3\.  
* **The Receiver:** While the message stack is being processed, the receiving processor must check for these incoming count messages 4\. Because ranks might not be synchronized on the same loop iteration, the receiver should **"aggressively drain the message stack"** to determine the total number of photon data messages ($N\_{mesg}$) it is expecting from all neighbors 4\.

### 4\. Preparation for Data Transfer

Once the handshake is complete and a rank knows its expected $N\_{mesg}$:

* **Post Receives:** The receiving rank posts exactly $N\_{mesg}$ calls to **MPI\_Irecv** 4\.  
* **Buffer Allocation:** For each of these MPI\_Irecv calls, a data buffer of size $N\_{max}$ must be pre-allocated to accommodate the incoming photon data 3\.

### 5\. Managing the Handshake in a Loop

To maintain a clean software design, this handshake should be part of a larger **asynchronous loop** 5\. If a processor has zero local rays to trace but is still expecting handshake messages (or subsequent data packets), it should use MPI\_Waitsome to wait for these counts to arrive 6\. If it still has local work, it uses MPI\_Testsome to check the handshake status without pausing computation 6\. The handshake is considered complete for an iteration only when all expected $N\_{mesg}$ values have been received and their corresponding data MPI\_Irecv calls have been posted 4, 6\.  
