#ifdef enablemultithread
#if (defined(__GNUC__))
#define TLS __thread
#elif (defined(_WIN32) || defined(_WIN64))
#define TLS __declspec(thread)
#endif
#else
#define TLS 
#endif

#ifdef enableatomic
#define ATOMICINT atomic_int
#else
#define ATOMICINT int 
#endif

extern TLS int commonAlloc1;
extern TLS int commonAlloc2;
extern TLS int **commonIP;
extern TLS int **commonJP;
