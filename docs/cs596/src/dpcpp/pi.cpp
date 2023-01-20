#include <CL/sycl.hpp>
#include <iostream>
#include <array>

using namespace cl::sycl;

#define NBIN 1000000  // # of bins for quadrature
#define NTRD 512      // # of threads

int main()
{
  float step = 1.0f/NBIN;
  std::array<float, NTRD> sum;
  for (int i=0; i<NTRD; ++i) sum[i] = 0.0f;

  queue q(gpu_selector{});

  std::cout << "Running on: " <<
  q.get_device().get_info<info::device::name>() << std::endl;

  range<1> sizeBuf{NTRD};

  {
    buffer<float, 1> sumBuf(sum.data(), sizeBuf);
    q.submit([&](handler &h){
      auto sumAccessor =
      sumBuf.get_access<access::mode::read_write>(h);
      h.parallel_for(sizeBuf, [=](id<1> tid) {
        for (int i=tid; i<NBIN; i+=NTRD) {
          float x = (i+0.5f)*step;
          sumAccessor[tid] += 4.0f/(1.0f+x*x);
        }
      });  // End parallel for
    });  // End queue submit
  }

  float pi=0.0f;
  for (int i=0; i<NTRD; i++)  // Inter-thread reduction
    pi += sum[i];
  pi *= step;  // Multiply bin width to complete integration

  std::cout << "Pi = " << pi << std::endl;

  return 0;
}

