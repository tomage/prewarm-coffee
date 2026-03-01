use leptos::*;

#[component]
pub fn AccordionGroup(
    title: &'static str,
    #[prop(default = false)]
    initially_open: bool,
    children: Children,
) -> impl IntoView {
    let (is_open, set_is_open) = create_signal(initially_open);

    view! {
        <div class="accordion-group" class:open=is_open>
            <button
                class="accordion-header"
                attr:aria-expanded=move || is_open.get().to_string()
                on:click=move |_| set_is_open.update(|v| *v = !*v)
            >
                <span class="accordion-title">{title}</span>
                <span class="accordion-arrow">{move || if is_open.get() { "\u{25B2}" } else { "\u{25BC}" }}</span>
            </button>
            <div class="accordion-body" class:collapsed=move || !is_open.get()>
                {children()}
            </div>
        </div>
    }
}
