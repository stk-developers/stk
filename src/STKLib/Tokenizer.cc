#include "Tokenizer.h"
#include <iostream>


namespace STK
{
  //****************************************************************************
  //****************************************************************************
  void
  Tokenizer::
  AddString(const char* pString)
  {
    // copy into string struct, which is more convenient
    std::string       aux_string(pString);
    std::string       aux_record;
    std::string::size_type    cur_pos = 0;
    std::string::size_type    old_pos = 0;

    // make sure we have enough space
    aux_record.reserve(aux_string.length());

    // find all of separators and make a list of tokens
    while(old_pos < std::string::npos) {
      // find the next separator
      cur_pos = aux_string.find_first_of(mSeparator, old_pos);
      std::cout << cur_pos << std::endl;

      // extract token
      aux_record.insert(0, pString, old_pos, cur_pos - old_pos);

      // insert to list
      this->push_back(aux_record);

      // update old position so that it points behind the separator
      old_pos = cur_pos < std::string::npos ? cur_pos + 1 : cur_pos;

      // we don't need the contents of the token
      aux_record.erase();
    }
  }


}; // namespace STK

