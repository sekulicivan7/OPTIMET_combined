// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#define CATCH_CONFIG_RUNNER
#define CATCH_CONFIG_CONSOLE_WIDTH 100

#include "Types.h"
#include "mpi/Session.h"
#include <catch.hpp>
#include <memory>
#include <random>
#include <regex>

std::unique_ptr<std::mt19937_64> mersenne(new std::mt19937_64(0));

int main(int argc, const char **argv) {
  Catch::Session session; // There must be exactly once instance
  optimet::mpi::init(argc, argv);

  // The following mess transforms the input arguments so that output files have different names
  // on different processors
  std::vector<std::string> arguments(argv, argv + argc);
  auto output_opt = std::find_if(arguments.begin(), arguments.end(), [](std::string const &arg) {
      if(arg == "-o" or arg == "--out")
      return true;
      auto const N = std::string("--out=").size();
      return arg.size() > N and arg.substr(0, N) == "--out=";
      });
  if(output_opt != arguments.end()) {
    if(*output_opt == "-o" or *output_opt == "--out")
      output_opt += 1;
    if(output_opt != arguments.end()) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if(rank > 0)
        *output_opt
          = std::regex_replace(*output_opt, std::regex("\\.xml"), std::to_string(rank) + ".xml");
    }
  }

  // transforms the modified arguments to a C-style array of pointers.
  std::vector<char const *> cargs(arguments.size());
  std::transform(arguments.begin(), arguments.end(), cargs.begin(),
      [](std::string const &c) { return c.c_str(); });

  int returnCode = session.applyCommandLine(argc, cargs.data());
  if(returnCode != 0) // Indicates a command line error
    return returnCode;
  mersenne.reset(new std::mt19937_64(session.configData().rngSeed));

  auto const result = session.run();
  optimet::mpi::finalize();
  return result;
}
