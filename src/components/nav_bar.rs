use leptos::prelude::*;
use leptos_router::components::A;
use leptos_router::hooks::use_location;

#[component]
pub fn NavBar() -> impl IntoView {
    let location = use_location();
    let pathname = move || location.pathname.get();

    let link_class = move |href: &'static str| {
        let current = pathname();
        if current == href || (href != "/" && current.starts_with(href)) {
            "nav-link active"
        } else {
            "nav-link"
        }
    };

    view! {
        <nav class="site-nav">
            <div class="site-nav-inner">
                <A href="/" attr:class="nav-brand">"Prewarm"</A>
                <div class="nav-links">
                    <A href="/" attr:class=move || link_class("/")>"Home"</A>
                    <A href="/theory" attr:class=move || link_class("/theory")>"Theory"</A>
                    <A href="/about" attr:class=move || link_class("/about")>"About"</A>
                </div>
            </div>
        </nav>
    }
}
