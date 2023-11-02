#include "console/chi_console.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#include <iostream>

//######################################################### Run Console loop
/** Executes the loop for the console.*/
void chi::Console::RunConsoleLoop(char*) const
{
  Chi::log.Log() << "Console loop started. "
                     << "Type \"exit\" to quit (or Ctl-C).";

  /** Wrapper to an MPI_Bcast call for a single integer
   * broadcast from location 0. */
  auto BroadcastSingleInteger = [](int* int_being_bcast)
  {
    MPI_Bcast(int_being_bcast,    //buffer
              1, MPI_INT,         //count + type
              0,                  //root
              Chi::mpi.comm);    //communicator
  };

  /** Wrapper to an MPI_Bcast call for an array of characters
   * broadcast from location 0. */
  auto HomeBroadcastStringAsRaw = [](std::string string_to_bcast,int length)
  {
    char* raw_string_to_bcast = string_to_bcast.data();
    MPI_Bcast(raw_string_to_bcast,           //buffer
              length, MPI_CHAR,              //count + type
              0,                             //root
              Chi::mpi.comm);               //communicator
  };

  /** Wrapper to an MPI_Bcast call for an array of characters
   * broadcast from location 0. This call is for non-home locations. */
  auto NonHomeBroadcastStringAsRaw = [](std::string& string_to_bcast,int length)
  {
    std::vector<char> raw_chars(length+1,'\0');
    MPI_Bcast(raw_chars.data(),              //buffer
              length, MPI_CHAR,              //count + type
              0,                             //root
              Chi::mpi.comm);               //communicator

    string_to_bcast = std::string(raw_chars.data());
  };

  /** Executes a string within the lua-console. */
  auto LuaDoString = [this](const std::string& the_string)
  {
    bool error = luaL_dostring(console_state_, the_string.c_str());
    if (error)
    {
      Chi::log.LogAll() << lua_tostring(console_state_, -1);
      lua_pop(console_state_, 1);
    }
  };

  auto ConsoleInputNumChars = [](const std::string& input)
  {
    int L = static_cast<int>(input.size());
    if (input == std::string("exit")) L = -1;

    return L;
  };

  const bool HOME = Chi::mpi.location_id == 0;

  while (not Chi::run_time::termination_posted_)
  {
    std::string console_input;

    if (HOME) std::cin >> console_input; //Home will be waiting here

    int console_input_len = ConsoleInputNumChars(console_input);

    BroadcastSingleInteger(&console_input_len); //Non-Home locs wait here

    if (console_input_len < 0) break;
    else
      if (HOME) HomeBroadcastStringAsRaw(console_input, console_input_len);
      else      NonHomeBroadcastStringAsRaw(console_input, console_input_len);

    try { LuaDoString(console_input); }
    catch(const Chi::RecoverableException& e)
    {
      Chi::log.LogAllError() << e.what();
    }
    catch(const std::exception& e)
    {
      Chi::log.LogAllError() << e.what();
      Chi::Exit(EXIT_FAILURE);
    }
  }//while not termination posted

  Chi::run_time::termination_posted_ = true;

  Chi::log.Log() << "Console loop stopped successfully.";
}
