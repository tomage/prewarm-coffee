/// Convert a time series from seconds to minutes: (t_s, value) → (t_min, value).
pub fn seconds_to_minutes(data: &[(f64, f64)]) -> Vec<(f64, f64)> {
    data.iter().map(|(t, v)| (t / 60.0, *v)).collect()
}

/// The leading part of a `(time, value)` series up to and including
/// `end_time`: the points a chart of that duration plots. Metrics computed on
/// it describe the run the user sees.
pub fn series_up_to(data: &[(f64, f64)], end_time: f64) -> Vec<(f64, f64)> {
    data.iter()
        .take_while(|(t, _)| *t <= end_time)
        .copied()
        .collect()
}

/// Bundled comparison metrics for cold-vs-preheated mug.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComparisonData {
    pub avg_temp_cold: f64,
    pub avg_temp_hot: f64,
    pub temp_diff: f64,
    /// `None` when the coffee never drops below the threshold within the series.
    pub crossing_time_cold: Option<f64>,
    /// `None` when the coffee never drops below the threshold within the series.
    pub crossing_time_hot: Option<f64>,
}

/// How much longer the preheated mug stays drinkable than the cold one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ExtraTime {
    /// Both crossings are known: the exact difference.
    Exact(f64),
    /// The preheated coffee never crosses within the run: at least this much longer.
    AtLeast(f64),
}

/// Extra drinkable time of the preheated mug over the cold one, in the same
/// time unit as the inputs. `duration` is the length of the simulated run.
///
/// Returns `None` when there is nothing meaningful to report: the cold coffee
/// never crossed (no baseline), or the preheated one crossed no later than it.
pub fn extra_time(cold: Option<f64>, hot: Option<f64>, duration: f64) -> Option<ExtraTime> {
    let cold = cold?;
    let extra = match hot {
        Some(h) => ExtraTime::Exact(h - cold),
        None => ExtraTime::AtLeast(duration - cold),
    };
    match extra {
        ExtraTime::Exact(d) | ExtraTime::AtLeast(d) if d <= 0.0 => None,
        e => Some(e),
    }
}

/// Compute all five comparison metrics from cold and hot time series.
///
/// Time-series data should be `(time, temp)` pairs; time units are preserved
/// as-is in the crossing times (caller decides whether minutes or seconds).
pub fn compute_comparison(
    cold_data: &[(f64, f64)],
    hot_data: &[(f64, f64)],
    threshold: f64,
) -> ComparisonData {
    let avg_temp_cold = average_while_above(cold_data, threshold);
    let avg_temp_hot = average_during_window(hot_data, cold_data, threshold);
    let crossing_time_cold = find_crossing_time_series(cold_data, threshold);
    let crossing_time_hot = find_crossing_time_series(hot_data, threshold);
    ComparisonData {
        avg_temp_cold,
        avg_temp_hot,
        temp_diff: avg_temp_hot - avg_temp_cold,
        crossing_time_cold,
        crossing_time_hot,
    }
}

/// Find when a generic time series crosses below a threshold (linearly interpolated).
///
/// Works on `(time, value)` tuples. Time units are preserved as-is (caller decides
/// whether they are seconds, minutes, etc.).
///
/// Returns `Some(time)` at the crossing, `Some(first_time)` if the series already
/// starts below the threshold, and `None` if it never drops below it (or is empty).
pub fn find_crossing_time_series(data: &[(f64, f64)], threshold: f64) -> Option<f64> {
    for i in 1..data.len() {
        let (t_prev, v_prev) = data[i - 1];
        let (t_curr, v_curr) = data[i];
        if v_prev >= threshold && v_curr < threshold {
            let frac = (threshold - v_prev) / (v_curr - v_prev);
            return Some(t_prev + frac * (t_curr - t_prev));
        }
    }
    match data.first() {
        Some(&(t, v)) if v < threshold => Some(t),
        _ => None,
    }
}

/// Average of values while they stay above `threshold`.
///
/// Takes the leading segment of `data` where values are >= threshold and returns
/// the mean. Returns 0.0 if the first value is already below threshold.
pub fn average_while_above(data: &[(f64, f64)], threshold: f64) -> f64 {
    let drinkable: Vec<_> = data.iter().take_while(|(_, v)| *v >= threshold).collect();
    if drinkable.is_empty() {
        return 0.0;
    }
    let sum: f64 = drinkable.iter().map(|(_, v)| v).sum();
    sum / drinkable.len() as f64
}

/// Average of `data` values over the window where `reference` stays above `threshold`.
///
/// Counts how many leading points in `reference` are >= threshold, then averages
/// the same number of leading points from `data`. Returns 0.0 if the reference
/// window is empty.
pub fn average_during_window(data: &[(f64, f64)], reference: &[(f64, f64)], threshold: f64) -> f64 {
    let ref_count = reference
        .iter()
        .take_while(|(_, v)| *v >= threshold)
        .count();
    if ref_count == 0 {
        return 0.0;
    }
    let sum: f64 = data.iter().take(ref_count).map(|(_, v)| v).sum();
    sum / ref_count as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- seconds_to_minutes ---

    #[test]
    fn test_seconds_to_minutes_basic() {
        let data = vec![(0.0, 95.0), (60.0, 90.0), (120.0, 85.0)];
        let result = seconds_to_minutes(&data);
        assert_eq!(result.len(), 3);
        assert!((result[0].0 - 0.0).abs() < 1e-10);
        assert!((result[1].0 - 1.0).abs() < 1e-10);
        assert!((result[2].0 - 2.0).abs() < 1e-10);
        // Values unchanged
        assert!((result[0].1 - 95.0).abs() < 1e-10);
        assert!((result[1].1 - 90.0).abs() < 1e-10);
    }

    #[test]
    fn test_seconds_to_minutes_empty() {
        let result = seconds_to_minutes(&[]);
        assert!(result.is_empty());
    }

    // --- find_crossing_time_series ---

    #[test]
    fn test_crossing_exact_boundary() {
        // Value hits threshold exactly at a data point
        let data = vec![(0.0, 80.0), (1.0, 60.0), (2.0, 50.0)];
        let t = find_crossing_time_series(&data, 60.0).unwrap();
        // 60.0 >= 60.0 at t=1, then 50 < 60 at t=2 → interpolation yields exactly t=1.0
        assert!((t - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_crossing_interpolation() {
        // Linear drop from 80 to 40 over 1 minute — crosses 60 at t=0.5
        let data = vec![(0.0, 80.0), (1.0, 40.0)];
        let t = find_crossing_time_series(&data, 60.0).unwrap();
        assert!((t - 0.5).abs() < 1e-10);
    }

    #[test]
    fn crossing_none_when_never_below() {
        let data = vec![(0.0, 90.0), (1.0, 80.0), (2.0, 70.0)];
        assert_eq!(find_crossing_time_series(&data, 60.0), None);
    }

    #[test]
    fn test_crossing_starts_below() {
        // Already below threshold at the first sample — crossed at the start
        let data = vec![(0.0, 50.0), (1.0, 40.0)];
        let t = find_crossing_time_series(&data, 60.0).unwrap();
        assert!(t.abs() < 1e-10);
    }

    #[test]
    fn test_crossing_empty_data() {
        let t = find_crossing_time_series(&[], 60.0);
        assert!(t.is_none());
    }

    #[test]
    fn test_crossing_single_point() {
        let data = vec![(5.0, 70.0)];
        assert_eq!(find_crossing_time_series(&data, 60.0), None);
    }

    // --- average_while_above ---

    #[test]
    fn test_avg_while_above_basic() {
        let data = vec![
            (0.0, 90.0),
            (1.0, 80.0),
            (2.0, 70.0),
            (3.0, 50.0),
            (4.0, 40.0),
        ];
        let avg = average_while_above(&data, 60.0);
        // Leading segment: 90, 80, 70 → mean = 80
        assert!((avg - 80.0).abs() < 1e-10);
    }

    #[test]
    fn test_avg_while_above_all_above() {
        let data = vec![(0.0, 90.0), (1.0, 80.0), (2.0, 70.0)];
        let avg = average_while_above(&data, 60.0);
        assert!((avg - 80.0).abs() < 1e-10);
    }

    #[test]
    fn test_avg_while_above_first_below() {
        let data = vec![(0.0, 50.0), (1.0, 80.0)];
        let avg = average_while_above(&data, 60.0);
        assert!((avg - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_avg_while_above_empty() {
        let avg = average_while_above(&[], 60.0);
        assert!((avg - 0.0).abs() < 1e-10);
    }

    // --- average_during_window ---

    #[test]
    fn test_avg_during_window_basic() {
        let reference = vec![(0.0, 90.0), (1.0, 80.0), (2.0, 70.0), (3.0, 50.0)];
        let data = vec![(0.0, 95.0), (1.0, 92.0), (2.0, 88.0), (3.0, 75.0)];
        let avg = average_during_window(&data, &reference, 60.0);
        // Reference has 3 leading points above 60 → average first 3 of data: (95+92+88)/3 = 91.666...
        assert!((avg - 91.666666).abs() < 0.001);
    }

    #[test]
    fn test_avg_during_window_ref_all_below() {
        let reference = vec![(0.0, 50.0), (1.0, 40.0)];
        let data = vec![(0.0, 95.0), (1.0, 90.0)];
        let avg = average_during_window(&data, &reference, 60.0);
        assert!((avg - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_avg_during_window_empty_reference() {
        let data = vec![(0.0, 95.0)];
        let avg = average_during_window(&data, &[], 60.0);
        assert!((avg - 0.0).abs() < 1e-10);
    }

    // --- compute_comparison ---

    #[test]
    fn test_compute_comparison_basic() {
        let cold = vec![(0.0, 90.0), (1.0, 75.0), (2.0, 55.0), (3.0, 40.0)];
        let hot = vec![(0.0, 90.0), (1.0, 80.0), (2.0, 65.0), (3.0, 50.0)];
        let result = compute_comparison(&cold, &hot, 60.0);

        // Cold crossing: between t=1 (75) and t=2 (55) → interpolated
        let cold_t = result.crossing_time_cold.unwrap();
        let hot_t = result.crossing_time_hot.unwrap();
        assert!(cold_t > 1.0 && cold_t < 2.0);
        // Hot crossing: between t=2 (65) and t=3 (50) → interpolated
        assert!(hot_t > 2.0 && hot_t < 3.0);
        // Hot avg should be >= cold avg (preheated mug keeps coffee warmer)
        assert!(result.temp_diff >= 0.0);
        assert!((result.temp_diff - (result.avg_temp_hot - result.avg_temp_cold)).abs() < 1e-10);
    }

    #[test]
    fn test_compute_comparison_identical_series() {
        let data = vec![(0.0, 90.0), (1.0, 70.0), (2.0, 50.0)];
        let result = compute_comparison(&data, &data, 60.0);
        assert!((result.temp_diff).abs() < 1e-10);
        assert!(
            (result.crossing_time_cold.unwrap() - result.crossing_time_hot.unwrap()).abs() < 1e-10
        );
    }

    #[test]
    fn comparison_carries_none() {
        let cold = vec![(0.0, 90.0), (1.0, 40.0)];
        let hot = vec![(0.0, 90.0), (1.0, 80.0)];
        let result = compute_comparison(&cold, &hot, 60.0);
        assert!(result.crossing_time_cold.is_some());
        assert_eq!(result.crossing_time_hot, None);
    }

    // --- extra_time ---

    #[test]
    fn extra_time_exact_when_both_cross() {
        assert_eq!(
            extra_time(Some(10.0), Some(14.0), 30.0),
            Some(ExtraTime::Exact(4.0))
        );
    }

    #[test]
    fn extra_time_lower_bound_when_hot_never_crosses() {
        assert_eq!(
            extra_time(Some(10.0), None, 30.0),
            Some(ExtraTime::AtLeast(20.0))
        );
    }

    #[test]
    fn extra_time_none_without_cold_baseline_or_gain() {
        assert_eq!(extra_time(None, None, 30.0), None);
        assert_eq!(extra_time(None, Some(5.0), 30.0), None);
        assert_eq!(extra_time(Some(10.0), Some(10.0), 30.0), None);
        assert_eq!(extra_time(Some(30.0), None, 30.0), None);
    }

    // --- series_up_to ---

    /// Linear per-minute series from `start` falling `slope` °C/min over
    /// `minutes` minutes.
    fn falling(start: f64, slope: f64, minutes: usize) -> Vec<(f64, f64)> {
        (0..=minutes)
            .map(|m| (m as f64, start - slope * m as f64))
            .collect()
    }

    #[test]
    fn series_up_to_keeps_points_within_duration() {
        let data = falling(90.0, 1.0, 120);
        let trimmed = series_up_to(&data, 30.0);
        assert_eq!(trimmed.len(), 31);
        assert_eq!(trimmed.last(), Some(&(30.0, 60.0)));
        assert!(series_up_to(&[], 30.0).is_empty());
    }

    #[test]
    fn metrics_on_trimmed_series_describe_the_displayed_run() {
        // Simulated for 120 min, displayed for 30. The preheated coffee
        // crosses at 117 min, after the displayed run ends.
        let duration = 30.0;
        let cold = series_up_to(&falling(100.0, 5.0, 120), duration); // crosses 50 at 10 min
        let hot = series_up_to(&falling(108.5, 0.5, 120), duration); // crosses 50 at 117 min
        let m = compute_comparison(&cold, &hot, 50.0);
        assert_eq!(m.crossing_time_cold, Some(10.0));
        assert_eq!(m.crossing_time_hot, None);
        assert_eq!(
            extra_time(m.crossing_time_cold, m.crossing_time_hot, duration),
            Some(ExtraTime::AtLeast(20.0))
        );
    }
}
