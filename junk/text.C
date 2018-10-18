#include <limits>
#include <iostream>

int main(int argc, char **argv)
{
  std::cout << "double max: " << std::numeric_limits<double>::max() << "\n"
            << "double min: " << std::numeric_limits<double>::min() << "\n";
  std::cout << "float max: " << std::numeric_limits<float>::max() << "\n"
            << "float min: " << std::numeric_limits<float>::min() << "\n";

  std::less<double> cmp;
  std::cout << "cmp(0,1) for std::less<double> " << cmp(0,1) << "\n";
  std::cout << "try this: "
            << (cmp(0,1) ?  std::numeric_limits<double>::max()
                         : -1.0 * std::numeric_limits<double>::min()) << "\n";
 
}
