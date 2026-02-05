use leptos::*;
use leptos_router::A;

#[component]
pub fn TheoryPage() -> impl IntoView {
    view! {
        <main class="container theory-page">
            <header>
                <h1>"Theory"</h1>
                <p class="tagline">"The physics behind the simulation"</p>
            </header>

            <nav class="back-nav">
                <A href="/">"< Back home"</A>
            </nav>

            <section class="theory-intro">
                <p>
                    "The " <A href="/">"simulator"</A> " models coffee cooling with two physical models, "
                    "each solved both numerically and analytically\u{2014}four methods in total. "
                    "The " <A href="/?tab=compare">"Compare tab"</A> " on the home page runs all four side by side."
                </p>
            </section>

            <section class="theory-section">
                <h2>"Lumped Capacitance"</h2>
                <p>
                    "The simpler model: treat the coffee and mug as two lumps of thermal mass "
                    "exchanging heat with each other and the environment via Newton\u{2019}s Law of Cooling."
                </p>
                <p>
                    "Solved two ways: a forward-Euler numerical integration and an "
                    "analytical closed-form solution. Good for intuition and fast to compute."
                </p>
                <p class="theory-link">
                    <A href="/theory/lumped">"See the full model and formulas >"</A>
                </p>
            </section>

            <section class="theory-section">
                <h2>"Heat Equation (PDE)"</h2>
                <p>
                    "A more realistic model: solve the radial heat equation in spherical coordinates "
                    "across three regions\u{2014}liquid, mug wall, and surrounding air. "
                    "Temperature varies continuously through each layer, not just between them."
                </p>
                <p>
                    "Solved two ways: a Crank\u{2013}Nicolson finite-difference scheme and an "
                    "analytical eigenfunction expansion. Captures spatial gradients\u{2014}the center "
                    "of the coffee stays hotter than the edges."
                </p>
                <p class="theory-link">
                    <A href="/theory/heat-equation">"See the full model and formulas >"</A>
                </p>
            </section>

            <section class="theory-section">
                <p>
                    "Want to know more about the project? Visit the "
                    <A href="/about">"about page"</A>
                    "."
                </p>
            </section>

            <nav class="back-nav bottom">
                <A href="/">"< Back home"</A>
            </nav>

            <footer>
                <p>
                    <a href="https://github.com/tomage/prewarm-coffee" target="_blank" rel="noopener">"Source code"</a>
                    " · AGPL-3.0 · "
                    "2025"
                </p>
            </footer>
        </main>
    }
}
