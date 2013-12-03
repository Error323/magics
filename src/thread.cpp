#include "thread.h"

static void *runThread(void *arg)
{
  return ((Thread *)arg)->run();
}

Thread::Thread():
  mTid(0),
  mRunning(0),
  mDetached(0) {}

Thread::~Thread()
{
  if (mRunning == 1 && mDetached == 0)
    pthread_detach(mTid);

  if (mRunning == 1)
    pthread_cancel(mTid);
}

int Thread::start()
{
  int result = pthread_create(&mTid, NULL, runThread, this);

  if (result == 0)
    mRunning = 1;

  return result;
}

int Thread::join()
{
  int result = -1;

  if (mRunning == 1)
  {
    result = pthread_join(mTid, NULL);

    if (result == 0)
      mDetached = 0;
  }

  return result;
}

int Thread::detach()
{
  int result = -1;

  if (mRunning == 1 && mDetached == 0)
  {
    result = pthread_detach(mTid);

    if (result == 0)
      mDetached = 1;
  }

  return result;
}

pthread_t Thread::self()
{
  return mTid;
}
