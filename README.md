# CP-scheduling
constraint-based scheduling code for MPTCP linux implementation version 0.89

This is patch code of Constraint-based Proactive Scheduling (CP-scheduling) for MPTCP linux implementation.
Currently, this patch version is incomplete so it may be suitable for research state. 
In current status, Note that this code may be only vailed for the 2 path case. (currently i will plan to conduct various tests for the CP-scheduling with more than 2 paths)

In order to employ the CP-scheduling, overwrite the patch code on the linux implementation as the following steps.

1. overwrite the mptcp_sched.c file to mptcp/net/mptcp/mptcp_sched.c file
2. overwirite the tcp.h file in the folder (mptcp_include_linux) to mptcp/include/linux/tcp.h file
3. update the modified the MPTCP linux kernel.
 
Then, you can use CP-scheduling when mptcp enables 'fullmech' state. This means default version is cannot applied after the patch is applied.

Bong-Hwan and Oh Jaiyong Lee, all rights reserved




