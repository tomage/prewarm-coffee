//! Browser adapter for [`Throttle`]: a signal that follows a memo at most once
//! per window, with a trailing update so it always settles on the latest value.

use crate::timing::{Throttle, ThrottleAction};
use leptos::leptos_dom::helpers::TimeoutHandle;
use leptos::*;
use std::cell::{Cell, RefCell};
use std::rc::Rc;
use std::time::Duration;

/// Returns a signal that takes `source`'s value on the leading edge of a burst
/// of changes, at most once per `window_ms` while the burst lasts, and once
/// more after it ends so it always rests on `source`'s latest value.
pub fn throttled<T: Clone + PartialEq + 'static>(source: Memo<T>, window_ms: f64) -> ReadSignal<T> {
    let (output, set_output) = create_signal(source.get_untracked());
    let throttle = Rc::new(RefCell::new(Throttle::new(window_ms)));
    // One timer at a time. `Throttle::on_timer` does not check the deadline,
    // so a leading run clears any armed timer rather than letting a stale one
    // fire a trailing run early.
    let timer: Rc<Cell<Option<TimeoutHandle>>> = Rc::new(Cell::new(None));

    let publish = move |value: T| {
        if output.with_untracked(|current| *current != value) {
            set_output.set(value);
        }
    };

    {
        let timer = timer.clone();
        create_effect(move |prev: Option<()>| {
            let value = source.get();
            if prev.is_none() {
                // Initial run: the signal was seeded from the source already.
                publish(value);
                return;
            }
            let now = js_sys::Date::now();
            let action = throttle.borrow_mut().on_input(now);
            match action {
                ThrottleAction::RunNow => {
                    if let Some(handle) = timer.take() {
                        handle.clear();
                    }
                    publish(value);
                }
                ThrottleAction::ScheduleTrailing { at_ms } => {
                    let throttle = throttle.clone();
                    let fired = timer.clone();
                    let delay = Duration::from_secs_f64((at_ms - now).max(0.0) / 1000.0);
                    let handle = set_timeout_with_handle(
                        move || {
                            fired.set(None);
                            if throttle.borrow_mut().on_timer(js_sys::Date::now()) {
                                publish(source.get_untracked());
                            }
                        },
                        delay,
                    )
                    .expect("set_timeout is available in the browser");
                    timer.set(Some(handle));
                }
                ThrottleAction::Wait => {}
            }
        });
    }

    on_cleanup(move || {
        if let Some(handle) = timer.take() {
            handle.clear();
        }
    });

    output
}
