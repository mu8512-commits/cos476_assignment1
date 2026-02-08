#include <stdio.h>
#include <thread>
#include <cstdlib>

#include "CycleTimer.h"

typedef struct {
    float x0, x1;
    float y0, y1;
    unsigned int width;
    unsigned int height;
    int maxIterations;
    int* output;
    int threadId;
    int numThreads;
} WorkerArgs;


extern void mandelbrotSerial(
    float x0, float y0, float x1, float y1,
    int width, int height,
    int startRow, int numRows,
    int maxIterations,
    int output[]);


void workerThreadStart_v1(WorkerArgs * const args) {

    double startTime = CycleTimer::currentSeconds();

    const int tid = args->threadId;
    const int nt  = args->numThreads;
    const int H   = (int)args->height;

    const int block = H / nt;          // floor
    int startRow = tid * block;
    int numRows  = block;

    // Last thread takes the remainder
    if (tid == nt - 1) {
        numRows = H - startRow;
    }

    // Optional safety: if more threads than rows, some threads do nothing
    if (numRows <= 0) return;

    mandelbrotSerial(args->x0, args->y0, args->x1, args->y1,
                     (int)args->width, (int)args->height,
                     startRow, numRows,
                     args->maxIterations,
                     args->output);
}

//
// workerThreadStart --
//
// Thread entrypoint.
void workerThreadStart_v2(WorkerArgs * const args) {

    double startTime = CycleTimer::currentSeconds();

    const int tid = args->threadId;
    const int nt  = args->numThreads;
    const int H   = (int)args->height;

    // Tune this: smaller => better balance, more overhead.
    // Start with 1 or 2 for Mandelbrot rows.
    const int CHUNK = 1;

    for (int startRow = tid * CHUNK; startRow < H; startRow += nt * CHUNK) {

        int numRows = CHUNK;
        if (startRow + numRows > H) numRows = H - startRow;
        if (numRows <= 0) break;

        mandelbrotSerial(args->x0, args->y0, args->x1, args->y1,
                         (int)args->width, (int)args->height,
                         startRow, numRows,
                         args->maxIterations,
                         args->output);
    }
}

//
// MandelbrotThread --
//
// Multi-threaded implementation of mandelbrot set image generation.
// Threads of execution are created by spawning std::threads.
void mandelbrotThread(
    int implThread,
    int numThreads,
    float x0, float y0, float x1, float y1,
    int width, int height,
    int maxIterations, int output[])
{
    static constexpr int MAX_THREADS = 64;

    if (numThreads > MAX_THREADS)
    {
        fprintf(stderr, "Error: Max allowed threads is %d\n", MAX_THREADS);
        exit(1);
    }

    // Creates thread objects that do not yet represent a thread.
    std::thread workers[MAX_THREADS];
    WorkerArgs args[MAX_THREADS];

    for (int i=0; i<numThreads; i++) {
      
        args[i].x0 = x0;
        args[i].y0 = y0;
        args[i].x1 = x1;
        args[i].y1 = y1;
        args[i].width = width;
        args[i].height = height;
        args[i].maxIterations = maxIterations;
        args[i].numThreads = numThreads;
        args[i].output = output;
      
        args[i].threadId = i;
    }

    // Spawn the worker threads. Note that only numThreads-1 std::threads
    // are created and the main application thread is used as a worker
    // as well.
     
    if(implThread == 0){
	    for (int i=1; i<numThreads; i++) {
		    workers[i] = std::thread(workerThreadStart_v1, &args[i]);
	    }

	    workerThreadStart_v1(&args[0]);
    }
    else if(implThread == 1){
	    for (int i=1; i<numThreads; i++) {
		    workers[i] = std::thread(workerThreadStart_v2, &args[i]);
	    }

	    workerThreadStart_v2(&args[0]);
    } else {
        // Avoid joining threads when implThread not selected.
        return;
    }

    // join worker threads
    for (int i=1; i<numThreads; i++) {
        workers[i].join();
    }
}

