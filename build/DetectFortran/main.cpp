#include "fcmangling.h"
extern "C" {
void OPTIMET_FC_GLOBAL(pdgemr2d, PDGEMR2D)(int *m, int *n, double *A,
            int *IA, int *JA, int *descA, double *B, int *IB, int *JB,
            int *descB, int *gcontext);
}
int main(int nargs, char **argv) {
  OPTIMET_FC_GLOBAL(pdgemr2d, PDGEMR2D)(nullptr, nullptr, nullptr,
            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
  return 0;
}
