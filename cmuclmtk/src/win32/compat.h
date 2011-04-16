#ifndef COMPAT_H
#define COMPAT_H

#ifdef _WIN32
#define popen(x, y) _popen((x), (y))
#define pclose(x) _pclose((x))
#define srandom(x) srand((x))
#define random() rand()
#define mkdtemp(x) rr_mkdtemp((x))


/* For mktemp and mkdir */
#include <io.h>
#include <direct.h>

#endif

#endif
