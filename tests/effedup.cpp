// vim: set sw=4 ts=4 sts=4 expandtab :

#include "effedup.hpp"

namespace EffedUp {

void prnFail(const std::string &message)
{
    std::cout
        << "### Check has FAILED!\n"
		<< "### " << message << std::endl;
}

} // namespace EffedUp
