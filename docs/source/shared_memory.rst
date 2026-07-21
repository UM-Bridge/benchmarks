.. _shared-memory:
=============================================
Client-Server Communication via Shared Memory
=============================================

UM-Bridge uses HTTP to communicate between a UQ client and a model server
by default. However, this is not the most performant option (in HPC) as it 
carries around 1 millisecond of overhead per request. We provide a quicker
method that uses shared memory (RAM) to transfer data.

Currently, only the C++ and Python implementations are supported; instructions 
are documented in the respective sections below. When enabled, UM-Bridge will 
first perform a data exchange test, and will fallback to HTTP if this fails.


Enabling Shared Memory Communication in C++
===========================================



Enabling Shared Memory Communication in Python
==============================================