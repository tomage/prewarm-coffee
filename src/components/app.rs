use crate::components::{ComparePanel, HeatEquationPanel, LumpedPanel};
use crate::storage::SharedSettings;
use leptos::*;
use leptos_router::{A, use_query_map};

#[component]
pub fn App() -> impl IntoView {
    let settings = SharedSettings::instance();
    let active_tab = settings.active_tab;

    // Reactively track ?tab= query parameter — works for both direct loads
    // and SPA navigation (the route view may be kept alive, so we can't rely
    // on the component body re-running).
    let query = use_query_map();
    create_effect(move |_| {
        let params = query.get();
        if let Some(tab) = params.get("tab") {
            if matches!(tab.as_str(), "lumped" | "heat_equation" | "compare") {
                active_tab.set(tab.to_string());
            }
        }
    });

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
                    "A ceramic mug at room temperature can pull several degrees from your coffee "
                    "in the first few minutes as it warms up."
                </p>
                <p>
                    "By preheating your mug with hot water first, the mug starts at a higher temperature "
                    "and steals less heat from your coffee. The result? Your coffee stays warmer, longer—"
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
                {move || {
                    let tab = active_tab.get();
                    match tab.as_str() {
                        "heat_equation" => view! { <HeatEquationPanel/> }.into_view(),
                        "compare" => view! { <ComparePanel/> }.into_view(),
                        _ => view! { <LumpedPanel/> }.into_view(),
                    }
                }}
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
                    " · No cookies · No personal data"
                </p>
                <p>
                    <a href="https://github.com/tomage/prewarm-coffee" target="_blank" rel="noopener">"Source code"</a>
                    " · AGPL-3.0 · "
                    "© 2025"
                </p>
            </footer>
        </main>
    }
}
