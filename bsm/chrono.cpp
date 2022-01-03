#include "chrono.h"

namespace bsm::chrono {
    using namespace bsm::chrono;
    frac_years time_between(datetime const &initial, datetime const &final) {
        auto dur = final - initial;
        return duration_cast<frac_years>(dur);
    }

}
