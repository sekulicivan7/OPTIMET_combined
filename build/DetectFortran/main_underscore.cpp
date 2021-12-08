#include "fcmangling.h"
extern "C" {
void OPTIMET_FC_GLOBAL_(blacs_setup, BLACS_SETUP)(int *nprow, int *npcol);
}
int main(int nargs, char **argv) {
  OPTIMET_FC_GLOBAL_(blacs_setup, BLACS_SETUP)(nullptr, nullptr);
  return 0;
}
