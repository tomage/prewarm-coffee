//! Unit-test support for code that runs effects. Effects run on the async
//! executor's next tick, not when a signal is set, so a test creates its
//! signals and effects under an owner and runs the executor before asserting
//! on what an effect did.

use any_spawner::Executor;
use leptos::prelude::Owner;

/// Install the test executor (once per process) and make a fresh owner the
/// current one on this thread. Keep it alive for the test: dropping it
/// disposes the test's signals and effects.
pub fn test_owner() -> Owner {
    // Ignore AlreadySet: tests share the process-wide executor.
    let _ = Executor::init_futures_executor();
    let owner = Owner::new();
    owner.set();
    owner
}

/// Run every effect that is due, until none is: the point where, in the
/// browser, the current task has ended and the queued effects have run.
pub fn tick() {
    Executor::poll_local();
}
