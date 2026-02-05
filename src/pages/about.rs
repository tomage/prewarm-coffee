use leptos::*;
use leptos_router::A;

#[component]
pub fn AboutPage() -> impl IntoView {
    view! {
        <main class="container about-page">
            <header>
                <h1>"Prewarm"</h1>
                <p class="tagline">"About this project"</p>
            </header>

            <nav class="back-nav">
                <A href="/">"< Back home"</A>
            </nav>

            <section class="about-section">
                <h2>"The Story"</h2>
                <p>
                    "Every time I make coffee, I warm up my coffee mug. It's perhaps somewhat "
                    "obvious that a hot mug keeps the coffee warmer longer, but I've always "
                    "wanted to know, for sure, how much it mattered, exactly. I could have "
                    "gone the empirical way and actually measured it a few times, and made "
                    "a little white-paper about it (I still might, come to think about it) "
                    "but, knowing that this can be solved both analytically and numerically, "
                    "I decided\u{2014}why not make a fun project out of it, and learn something "
                    "along the way?"
                </p>
            </section>

            <section class="about-section">
                <h2>"The Build"</h2>
                <p>
                    "The last time I modeled a real physical system was a heat transfer problem "
                    "about boiling an egg, back in university. I'd been wanting to do something "
                    "like that again ever since. I was also curious about WebAssembly\u{2014}what "
                    "it feels like to build with, how it performs in the browser. Coffee cooling "
                    "turned out to be the perfect excuse to do both at the same time."
                </p>
                <p>
                    "Everything runs client-side in your browser using Rust compiled to WASM. "
                    "No servers, no cookies, no personal data\u{2014}just physics "
                    "(and some privacy-friendly analytics)."
                </p>
                <p>
                    "Curious about the math? Check out the "
                    <A href="/theory">"theory and models"</A>
                    " behind the simulation, or head back to the "
                    <A href="/">"interactive simulator"</A>
                    " to experiment yourself."
                </p>
            </section>

            <section class="about-section">
                <h2>"Version History"</h2>
                <dl class="version-history">
                    <dt>"v2 \u{2014} Current " <span class="version-date">"February 8, 2026"</span></dt>
                    <dd>
                        <ul>
                            <li>"Added a heat-equation model: 3-region PDE across liquid, mug wall, and air"</li>
                            <li>"Improved lumped model with analytical solvers, evaporation, radiation, and lid-opening controls"</li>
                            <li>"Interactive profile charts, comparison metrics, and theory pages"</li>
                        </ul>
                    </dd>
                    <dt>"v1 \u{2014} Initial release " <span class="version-date">"February 4, 2026"</span></dt>
                    <dd>
                        <ul>
                            <li>"Single lumped-capacitance numerical solver (forward Euler)"</li>
                            <li>"Basic cooling chart comparing preheated vs cold mug"</li>
                            <li>"Material selection, volume, and wall thickness sliders"</li>
                        </ul>
                    </dd>
                </dl>
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
