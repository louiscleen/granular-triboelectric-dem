#ifndef PROGRESS_TRACKER_HPP
#define PROGRESS_TRACKER_HPP

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

class ProgressTracker {
  public:
    static constexpr std::chrono::steady_clock::duration default_report_interval =
        std::chrono::minutes(10);

    // Note : pas besoin de template pour report_interval_ car <chrono> gère déjà les conversions
    // entre durées
    ProgressTracker(int tasks_total, int n_workers, std::ostream &os,
                    std::chrono::steady_clock::duration report_interval = default_report_interval,
                    double alpha = 0.001)
        : tasks_total_(tasks_total), n_workers_(n_workers), os_(os),
          report_interval_(report_interval), alpha_(alpha), task_duration_per_worker_(n_workers),
          sizes_(n_workers, 0) {
        print_report();
    }

    void update(int worker, std::chrono::steady_clock::duration task_duration) {
        task_duration_per_worker_[worker].push_back(task_duration);
        sizes_[worker]++;
        if (sizes_[worker] == min_size_ + 1) {
            min_size_ = *std::min_element(sizes_.begin(), sizes_.end());
        }

        tasks_completed_++;
        if (task_duration > max_task_duration_) {
            max_task_duration_ = task_duration;
        }
        compute_EMA(task_duration);
        if (should_print_report()) {
            print_report();
            time_at_last_report_ = std::chrono::steady_clock::now();
            tasks_completed_at_last_report_ = tasks_completed_;
        }
    }

  private:
    using steady_clock = std::chrono::steady_clock;
    using duration = steady_clock::duration;

    // ---- Configuration ----
    const int tasks_total_;
    const int n_workers_;
    const steady_clock::time_point start_time_ = steady_clock::now();
    std::ostream &os_;
    const duration report_interval_;
    const double alpha_;
    const double warmup_threshold_ =
        std::max(2.0 * n_workers_, 0.1 * tasks_total_); // Nombre de tâches minimum avant que l'on
                                                        // fournisse un ETA (car pas fiable avant)

    // ---- State ----
    steady_clock::time_point time_at_last_report_ = steady_clock::now();
    int tasks_completed_ = 0;
    int tasks_completed_at_last_report_ = 0;
    duration max_task_duration_ = duration::zero();
    mutable double ema_ = 0.0;
    std::vector<std::vector<duration>> task_duration_per_worker_;
    std::vector<std::size_t> sizes_;
    std::size_t min_size_ = 0;
    mutable std::size_t old_min_size_ = 0;

    void compute_EMA(duration task_duration) const {
        double task_duration_sec = std::chrono::duration<double>(task_duration).count();
        if (task_duration_sec <= 0.0)
            return;

        if (min_size_ == 0)
            return;

        if (old_min_size_ == min_size_ && tasks_total_ - tasks_completed_ > n_workers_)
            return;

        old_min_size_ = min_size_;

        double active_workers;
        if (tasks_total_ - tasks_completed_ < n_workers_)
            active_workers = tasks_total_ - tasks_completed_;
        else
            active_workers = n_workers_;

        double mean_duration_sec = 0.0;
        for (int i = 0; i < n_workers_; ++i) {
            mean_duration_sec +=
                std::chrono::duration<double>(task_duration_per_worker_[i][min_size_ - 1]).count();
        };
        mean_duration_sec /= n_workers_;

        double throughput = (1.0 / mean_duration_sec) * active_workers; // tasks per second
        if (ema_ == 0.0) {
            ema_ = throughput;
        } else {
            ema_ = alpha_ * throughput + (1.0 - alpha_) * ema_;
        }
    }

    std::optional<duration> compute_ETA() const {
        if (min_size_ == 0) // tasks_completed_ < warmup_threshold_ ||
            return std::nullopt;

        int tasks_remaining = tasks_total_ - tasks_completed_;
        double eta = tasks_remaining / ema_;
        return std::chrono::duration_cast<duration>(std::chrono::duration<double>{eta});
    }

    bool should_print_report() const {
        if (tasks_completed_ == tasks_total_) {
            return true;
        }
        return (steady_clock::now() - time_at_last_report_ >= report_interval_ &&
                tasks_completed_ > tasks_completed_at_last_report_);
    }

    void print_report() const {
        int progress_percent = static_cast<int>(100.0 * tasks_completed_ / tasks_total_);
        const int tasks_active = std::min(n_workers_, tasks_total_ - tasks_completed_);
        const auto elapsed = steady_clock::now() - start_time_;
        const auto eta = compute_ETA();

        os_ << "[" << std::setw(3) << progress_percent << "%] " << std::setw(5) << tasks_completed_
            << "/" << tasks_total_ << " | active: " << std::setw(4) << tasks_active
            << " | Elapsed = " << format_hhmmss_ms(elapsed)
            << " | Throughput (EMA) = " << std::fixed << std::setprecision(4) << (ema_)
            << " tasks/s"
            << " | ETA = ";
        if (eta) {
            os_ << format_hhmmss_ms(*eta) << "\n";

        } else {
            os_ << "estimating...\n";
        }
        os_.flush();
    }

    static std::string format_hhmmss(duration d) {

        std::chrono::hh_mm_ss hms(d);

        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(2) << hms.hours().count() << ":" << std::setw(2)
            << hms.minutes().count() << ":" << std::setw(2) << hms.seconds().count();
        return oss.str();
    }

    static std::string format_hhmmss_ms(duration d) {
        auto d_ms = std::chrono::duration_cast<std::chrono::milliseconds>(d);

        std::chrono::hh_mm_ss hms(d_ms);

        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(2) << hms.hours().count() << ":" << std::setw(2)
            << hms.minutes().count() << ":" << std::setw(2) << hms.seconds().count() << "."
            << std::setw(3) << hms.subseconds().count();

        return oss.str();
    }
};

#endif // PROGRESS_TRACKER_HPP

// DESTRUCTEUR ???
