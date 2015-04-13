#include <ctime>

int GetTime(int start = 0)
{
  static time_t _SAVETIME_;
  if (start) {
    _SAVETIME_=time(0);
    return 0;
  }
  else {
    return (time(0)-_SAVETIME_);
  }

}
