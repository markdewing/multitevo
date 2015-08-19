
from mpi4py import MPI

comm = None

MyRank = 0
NRanks = 1

def initParallel():
    global MyRank, NRanks, comm
    comm = MPI.COMM_WORLD
    MyRank = comm.Get_rank()
    NRanks = comm.Get_size()

def printRank():
    if MyRank == 0:
        return True
    return False

def sendReceiveParallel(sendBuf, dest, source):
    if not comm:
        return sendBuf
    return comm.sendrecv(sendBuf, dest, source=source)


def addParallel(sendBuf, recvBuf):
    global comm
    if not comm:
        recvBuf[:] = sendBuf[:]
        return
    comm.Allreduce(sendBuf, recvBuf)

