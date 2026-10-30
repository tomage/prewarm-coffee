//! Pure throttle and debounce state machines.
//!
//! No clocks or timers live here: callers pass `now_ms` in and schedule the
//! returned deadlines themselves, so the logic is testable without a browser.

/// Minimum spacing between expensive solves while a slider is dragged.
/// 100 ms keeps the chart live during a drag while capping solves at ~10/s.
pub const SOLVE_THROTTLE_MS: f64 = 100.0;

/// Quiet period before settings are persisted. 250 ms means a drag results in
/// one localStorage write after it settles, not one per input event.
pub const SAVE_DEBOUNCE_MS: f64 = 250.0;

/// Quiet period before a page's scroll position is saved. 100 ms writes once
/// per scroll gesture instead of once per scroll event, while a position held
/// for a moment is saved before the user can plausibly navigate away (leaving
/// the page also flushes the pending save).
pub const SCROLL_SAVE_DEBOUNCE_MS: f64 = 100.0;

/// What the caller should do after a throttled input.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThrottleAction {
    /// Run the work immediately (leading edge).
    RunNow,
    /// Arm a timer for `at_ms`; call [`Throttle::on_timer`] when it fires.
    ScheduleTrailing { at_ms: f64 },
    /// A trailing run is already scheduled; do nothing.
    Wait,
}

/// Leading + trailing throttle: at most one run per window, and the last
/// input of a burst is always followed by a run.
#[derive(Debug, Clone)]
pub struct Throttle {
    window_ms: f64,
    last_run_ms: Option<f64>,
    trailing_pending: bool,
}

impl Throttle {
    pub fn new(window_ms: f64) -> Self {
        Self {
            window_ms,
            last_run_ms: None,
            trailing_pending: false,
        }
    }

    pub fn on_input(&mut self, now_ms: f64) -> ThrottleAction {
        match self.last_run_ms {
            Some(last) if now_ms - last < self.window_ms => {
                if self.trailing_pending {
                    ThrottleAction::Wait
                } else {
                    self.trailing_pending = true;
                    ThrottleAction::ScheduleTrailing {
                        at_ms: last + self.window_ms,
                    }
                }
            }
            _ => {
                self.last_run_ms = Some(now_ms);
                self.trailing_pending = false;
                ThrottleAction::RunNow
            }
        }
    }

    /// Returns true if the trailing run should happen now; restarts the window.
    pub fn on_timer(&mut self, now_ms: f64) -> bool {
        if !self.trailing_pending {
            return false;
        }
        self.trailing_pending = false;
        self.last_run_ms = Some(now_ms);
        true
    }
}

/// Trailing debounce: fires once, `delay_ms` after the last input.
#[derive(Debug, Clone)]
pub struct Debounce {
    delay_ms: f64,
    deadline_ms: Option<f64>,
}

impl Debounce {
    pub fn new(delay_ms: f64) -> Self {
        Self {
            delay_ms,
            deadline_ms: None,
        }
    }

    /// Records an input and returns the (pushed) deadline to arm a timer for.
    pub fn on_input(&mut self, now_ms: f64) -> f64 {
        let deadline = now_ms + self.delay_ms;
        self.deadline_ms = Some(deadline);
        deadline
    }

    /// True only if an input is pending and its deadline has been reached.
    pub fn on_timer(&mut self, now_ms: f64) -> bool {
        match self.deadline_ms {
            Some(deadline) if now_ms >= deadline => {
                self.deadline_ms = None;
                true
            }
            _ => false,
        }
    }

    /// The deadline of the pending input, if any. A timer that fired before
    /// its deadline (timers may run early) re-arms for the time remaining.
    pub fn deadline(&self) -> Option<f64> {
        self.deadline_ms
    }

    /// Clears and reports a pending input, for a flush on `pagehide`.
    pub fn take_pending(&mut self) -> bool {
        self.deadline_ms.take().is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_input_runs_now() {
        let mut t = Throttle::new(100.0);
        assert_eq!(t.on_input(0.0), ThrottleAction::RunNow);
    }

    #[test]
    fn inputs_within_window_schedule_one_trailing() {
        let mut t = Throttle::new(100.0);
        assert_eq!(t.on_input(0.0), ThrottleAction::RunNow);
        assert_eq!(
            t.on_input(10.0),
            ThrottleAction::ScheduleTrailing { at_ms: 100.0 }
        );
        assert_eq!(t.on_input(50.0), ThrottleAction::Wait);
    }

    #[test]
    fn trailing_timer_runs_and_restarts_window() {
        let mut t = Throttle::new(100.0);
        t.on_input(0.0);
        t.on_input(10.0);
        assert!(t.on_timer(100.0));
        assert_eq!(
            t.on_input(110.0),
            ThrottleAction::ScheduleTrailing { at_ms: 200.0 }
        );
    }

    #[test]
    fn input_after_quiet_window_runs_now() {
        let mut t = Throttle::new(100.0);
        t.on_input(0.0);
        assert_eq!(t.on_input(250.0), ThrottleAction::RunNow);
    }

    #[test]
    fn no_trailing_when_no_input_after_leading() {
        let mut t = Throttle::new(100.0);
        t.on_input(0.0);
        assert!(!t.on_timer(100.0));
    }

    #[test]
    fn debounce_pushes_deadline() {
        let mut d = Debounce::new(250.0);
        assert_eq!(d.on_input(0.0), 250.0);
        assert_eq!(d.on_input(100.0), 350.0);
        assert!(!d.on_timer(250.0));
        assert!(d.on_timer(350.0));
    }

    #[test]
    fn debounce_fires_once_after_quiet() {
        let mut d = Debounce::new(250.0);
        assert!(!d.on_timer(0.0));
        d.on_input(0.0);
        assert!(d.on_timer(250.0));
        assert!(!d.on_timer(300.0));
    }

    #[test]
    fn debounce_early_timer_keeps_deadline() {
        let mut d = Debounce::new(250.0);
        assert_eq!(d.deadline(), None);
        d.on_input(0.0);
        assert!(!d.on_timer(249.5));
        assert_eq!(d.deadline(), Some(250.0));
        assert!(d.on_timer(250.0));
        assert_eq!(d.deadline(), None);
    }

    #[test]
    fn debounce_take_pending_flushes() {
        let mut d = Debounce::new(250.0);
        assert!(!d.take_pending());
        d.on_input(0.0);
        assert!(d.take_pending());
        assert!(!d.take_pending());
        assert!(!d.on_timer(500.0));
    }
}
