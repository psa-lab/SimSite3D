/******************************************************************************
 * Copyright (c) 2010, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/
#include <IK_tests.H>

using namespace ASCbase;

int main(const int argc, const char **argv)
{
  std::cout << "\n" << argv[0] << " (" << PACKAGE_NAME << ") " 
            << PACKAGE_VERSION << "\n\n";

  // Do not return -1 here since system, fork, etc return -1 on failure, and we
  // wish to distinguish between system and program failure
  SearchParameters my_params(argc, argv);
  BaseParameters::status_t status = my_params.status();
  if(status == BaseParameters::DISPLAY_HELP_ONLY) return 0;
  else if(status != BaseParameters::READY){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize parameters\n";
    return 1;
  }
  IK_tests my_search(&my_params);
  if(my_search.fail()){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize search\n";
    return 1;
  }

  my_search.test_moving_stuff();
}
