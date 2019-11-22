#include <sys/resource.h>
#include "getmem.h"

/* メモリ使用量取得(cから呼び出し) */
int getmem()
{
  int mem;
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  mem = ru.ru_maxrss;

  return mem;
}
