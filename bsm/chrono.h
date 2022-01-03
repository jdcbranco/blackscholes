#ifndef BSM_CHRONO_H
#define BSM_CHRONO_H

#include <concepts>
#include <chrono>
#include <vector>
#include <iostream>
#include <ratio>

namespace bsm::chrono {
    using namespace std::chrono;
    using namespace std::chrono_literals;

    using frac_years = std::chrono::duration<long double, std::ratio<31556952>>;
    using frac_seconds = std::chrono::duration<long double>;
    using frac_nanoseconds	= std::chrono::duration<long double, std::nano>;
    using frac_time = std::chrono::time_point<std::chrono::system_clock, frac_nanoseconds>;

    constexpr frac_years operator""_y(long double __y) noexcept {
        return frac_years{__y};
    }

    constexpr frac_years operator""_years(long double __y) noexcept {
        return frac_years{__y};
    }

    struct datetime2 {
        std::chrono::time_point<std::chrono::system_clock> instant;
        inline datetime2(std::chrono::year_month_day const& date, std::chrono::seconds time): instant{std::chrono::sys_days{date}+time} {
        }
        inline datetime2(std::chrono::year_month_day const& date): datetime2{std::chrono::sys_days{date},0s} {
        }
        inline datetime2(std::chrono::time_point<std::chrono::system_clock> instant): instant{instant} {
        }
        inline datetime2(frac_time instant): instant{std::chrono::duration_cast<seconds>(instant.time_since_epoch())} {
        }

        inline datetime2(datetime2 const&) = default;
        inline datetime2(datetime2 &&) noexcept = default;
        inline static datetime2 now() noexcept {
            return datetime2{std::chrono::system_clock::now()};
        }
        datetime2 operator+(std::chrono::seconds time) {
            datetime2 copy{*this};
            copy.instant += time;
            return copy;
        }
        datetime2 operator+(frac_seconds time) {
            datetime2 copy{*this};
            copy.instant += duration_cast<std::chrono::seconds>(time);
            return copy;
        }

        auto operator-(datetime2 const& other) const {
            return this->instant - other.instant;
        }

        bool operator==(datetime2 const& other) const {
            return this->instant == other.instant;
        }
    };

    using datetime = datetime2;

    frac_years time_between(datetime const& initial, datetime const& final);

//    auto time_between(datetime const& initial, datetime const& final) {
//        auto dur = final - initial;
//        return duration_cast<frac_years>(dur);
//    }

//    class calendar {
//        std::vector<std::chrono::year_month_day> holidays;
//    public:
//        int business_days_between() const {
//            return 0;
//        }
//        bool isBusinessDay(std::chrono::year_month_day ymd)
//        {
//            std::chrono::sys_days sd = ymd;
//            std::chrono::weekday wd = sd;
//            if (wd == std::chrono::Saturday || wd == std::chrono::Sunday)  // weekend
//                return false;
//            return true;
//        }
//    };
}


#endif //BSM_CHRONO_H
