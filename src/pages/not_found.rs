use leptos::*;
use leptos_router::A;

#[component]
pub fn NotFoundPage() -> impl IntoView {
    view! {
        <main class="container about-page">
            <header>
                <h1>"404"</h1>
                <p class="tagline">"Page not found"</p>
            </header>

            <nav class="back-nav">
                <A href="/">"< Back home"</A>
            </nav>
        </main>
    }
}
