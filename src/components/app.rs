use crate::components::overlay_history::forget_overlay_entry;
use crate::components::{use_media_query, SimTab, SimulatorOverlay, MOBILE_QUERY};
use crate::storage::SharedSettings;
use leptos::*;
use leptos_router::{use_query_map, A};

#[component]
pub fn App() -> impl IntoView {
    let settings = SharedSettings::instance();
    let active_tab = settings.active_tab;

    let (show_overlay, set_show_overlay) = create_signal(false);
    let open_sim_button = create_node_ref::<html::Button>();

    // Exactly one simulator is mounted: the inline one, or on small screens
    // the entry card and its overlay. Leaving the small layout closes the
    // overlay, so coming back shows the entry card rather than a stale overlay.
    // It leaves the overlay's history entry in place: stepping back over it
    // here would race the resize, so the next Back pops it instead.
    let is_mobile = use_media_query(MOBILE_QUERY);
    create_effect(move |_| {
        if !is_mobile.get() {
            set_show_overlay.set(false);
        }
    });

    // Leaving this page leaves any overlay entry below the new page's; Back
    // onto it is then an ordinary navigation (see `overlay_history`).
    on_cleanup(forget_overlay_entry);

    // A ?tab= link selects that tab. Set only when it differs: a set notifies
    // even when unchanged, and the save effect would then write at startup.
    let query = use_query_map();
    create_effect(move |_| {
        let params = query.get();
        if let Some(tab) = params.get("tab") {
            if matches!(tab.as_str(), "lumped" | "heat_equation" | "compare")
                && active_tab.with_untracked(|current| current != tab)
            {
                active_tab.set(tab.to_string());
            }
        }
    });

    let render_active_panel = move || SimTab::from_key(&active_tab.get()).panel(true);

    let mobile_entry = move || {
        view! {
            <div class="mobile-sim-entry">
                <p>"Adjust parameters and watch how preheating affects cooling in real time."</p>
                <button
                    class="open-sim-btn"
                    node_ref=open_sim_button
                    on:click=move |_| set_show_overlay.set(true)
                >"Open Simulator"</button>
            </div>
        }
    };

    let desktop_sim = move || {
        view! {
            <div class="desktop-sim">
                <div class="solver-toggle solver-toggle--3">
                    <button
                        class:active=move || active_tab.get() == "lumped"
                        on:click=move |_| active_tab.set("lumped".to_string())
                    >"Lumped Capacitance"</button>
                    <button
                        class:active=move || active_tab.get() == "heat_equation"
                        on:click=move |_| active_tab.set("heat_equation".to_string())
                    >"Heat Equation"</button>
                    <button
                        class:active=move || active_tab.get() == "compare"
                        on:click=move |_| active_tab.set("compare".to_string())
                    >"Compare"</button>
                </div>
                {render_active_panel}
            </div>
        }
    };

    view! {
        <main class="container">
            <header>
                <h1>"Prewarm"</h1>
                <p class="tagline">"Why preheating your mug keeps coffee warmer"</p>
            </header>

            <section class="intro">
                <h2>"The Science"</h2>
                <p>
                    "When you pour hot coffee into a cold mug, the mug absorbs heat from your coffee. "
                    "A cold 180\u{00A0}ml ceramic mug takes about 13\u{00A0}\u{00B0}C out of your coffee "
                    "as the two even out; a mug preheated to 80\u{00A0}\u{00B0}C takes only about 2\u{00A0}\u{00B0}C."
                </p>
                <p>
                    "By preheating your mug with hot water first, the mug starts at a higher temperature "
                    "and steals less heat from your coffee. The result? Your coffee stays warmer, longer\u{2014}"
                    "and you get more time in the ideal drinking range."
                </p>
            </section>

            <section class="interactive">
                <h2>"Try It Yourself"</h2>
                <p class="wip-note">
                    "This is a work in progress. The models use simplified physics (spherical geometry, "
                    "uniform properties) and may not accurately reflect real-world mug cooling. "
                    "See "<A href="/theory">"Theory & Models"</A>" for assumptions and limitations."
                </p>

                {move || if is_mobile.get() {
                    mobile_entry().into_view()
                } else {
                    desktop_sim().into_view()
                }}
            </section>

            // Numbers checked against both models by physics/claims_tests.rs;
            // change them together.
            <section class="conclusion">
                <h2>"What It Adds Up To"</h2>
                <p>
                    "For a typical ceramic mug, preheating gives you about "
                    <strong>"3\u{00BD} more minutes"</strong>
                    " above drinking temperature, and your coffee stays "
                    <strong>"5\u{2013}8\u{00A0}\u{00B0}C warmer"</strong>
                    " on average while you drink it. For a 180\u{00A0}ml cup that\u{2019}s about a third more time, "
                    "and the extra minutes barely depend on mug size: a big mug\u{2019}s head start is smaller, "
                    "but it fades more slowly."
                </p>
                <p>
                    "The mug\u{2019}s walls decide how much it matters. A thick-walled mug (8\u{00A0}mm) gains "
                    "about 6 minutes; a thin one (2\u{00A0}mm) only about one. Preheating helps only as much "
                    "as the mug can soak up heat."
                </p>
                <p class="conclusion-note">
                    "These numbers come from the two models on this page, which agree closely. "
                    "They haven\u{2019}t been checked against a measured cup yet."
                </p>
            </section>

            <section class="learn-more">
                <h2>"Learn More"</h2>
                <div class="learn-more-links">
                    <A href="/theory" class="learn-more-card">
                        <h3>"Theory & Models"</h3>
                        <p>"Explore the physics behind the simulation\u{2014}from simple lumped capacitance to analytical eigenfunction solutions."</p>
                    </A>
                    <A href="/about" class="learn-more-card">
                        <h3>"About"</h3>
                        <p>"The story behind this project and why we built it."</p>
                    </A>
                </div>
            </section>

            <footer>
                <p class="disclaimer">
                    "This simulation uses simplified physics models and approximations. "
                    "Real-world results may vary based on factors not modeled here (evaporation, lid usage, mug shape, etc.). "
                    <a href="https://github.com/tomage/prewarm-coffee/issues" target="_blank" rel="noopener">"Report inaccuracies"</a>
                </p>
                <p class="analytics-note">
                    "Analytics by "
                    <a href="https://umami.is" target="_blank" rel="noopener">"Umami"</a>
                    " \u{00B7} No cookies \u{00B7} No personal data"
                </p>
                <p>
                    <a href="https://github.com/tomage/prewarm-coffee" target="_blank" rel="noopener">"Source code"</a>
                    " \u{00B7} AGPL-3.0 \u{00B7} "
                    "\u{00A9} 2025"
                </p>
            </footer>
        </main>

        <Show when=move || is_mobile.get()>
            <SimulatorOverlay
                show=show_overlay
                set_show=set_show_overlay
                active_tab=active_tab
                return_focus=open_sim_button
            />
        </Show>
    }
}
