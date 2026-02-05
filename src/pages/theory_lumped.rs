use crate::components::LumpedPanel;
use leptos::*;
use leptos_router::A;

#[component]
pub fn TheoryLumpedPage() -> impl IntoView {
    view! {
        <main class="container details-page">
            <header>
                <h1>"Theory: Lumped Capacitance"</h1>
                <p class="tagline">"How the simulation works"</p>
            </header>

            <nav class="back-nav">
                <A href="/theory">"< Back to Theory"</A>
            </nav>

            <p class="wip-note">
                "This page is a work in progress. There may be discrepancies between "
                "the explanations below and what the simulation actually computes, as well "
                "as outright errors. If you spot something, please "
                <a href="https://github.com/tomage/prewarm-coffee/issues" target="_blank" rel="noopener">"open an issue"</a>"."
            </p>

            <section class="physics-detail">
                <h2>"Newton's Law of Cooling"</h2>
                <p>
                    "The foundation of our simulation is Newton's Law of Cooling, which states that "
                    "the rate of heat loss of a body is proportional to the temperature difference "
                    "between the body and its surroundings."
                </p>
                <div class="formula">
                    "\\[\\frac{dT}{dt} = -k(T - T_{\\text{env}})\\]"
                </div>
                <p class="formula-explanation">
                    "Where "
                    "\\(T\\)" " is the object's temperature, "
                    "\\(T_{\\text{env}}\\)" " is the environment temperature, "
                    "and " "\\(k\\)" " is a cooling constant that depends on the object's properties."
                </p>
            </section>

            <section class="physics-detail">
                <h2>"Heat Transfer Model"</h2>
                <p>
                    "Our simulation uses a lumped capacitance model, tracking five heat flows:"
                </p>

                <div class="heat-flow">
                    <h3>"1. Coffee to Mug (Conduction)"</h3>
                    <div class="formula">
                        "\\[Q = k_{\\text{eff}} \\cdot (T_{\\text{coffee}} - T_{\\text{mug}})\\]"
                    </div>
                    <p>
                        "Heat transfers from the hot coffee to the mug walls through direct contact. "
                        "The effective coefficient " "\\(k_{\\text{eff}} = 2\\,k_{\\text{material}}\\)" " is a simplified "
                        "model for coffee-to-mug heat transfer, where " "\\(k_{\\text{material}}\\)" " is the mug "
                        "material's thermal conductivity in W/(m\u{00B7}K)."
                    </p>
                </div>

                <div class="heat-flow">
                    <h3>"2. Coffee to Environment (Convection)"</h3>
                    <div class="formula">
                        "\\[Q = f_{\\text{open}} \\cdot h \\cdot A_{\\text{top}} \\cdot (T_{\\text{coffee}} - T_{\\text{room}})\\]"
                    </div>
                    <p>
                        "Heat escapes from the exposed coffee surface (top of the mug) through "
                        "natural convection. " "\\(h\\)" " is the convective heat transfer coefficient "
                        "(~20 W/m"<sup>"2"</sup>"K for natural convection with buoyancy), "
                        "\\(A_{\\text{top}}\\)" " is the top surface area, and "
                        "\\(f_{\\text{open}}\\)" " is the lid opening fraction (0 = fully covered, 1 = fully open)."
                    </p>
                </div>

                <div class="heat-flow">
                    <h3>"3. Mug to Environment (Convection)"</h3>
                    <div class="formula">
                        "\\[Q = h \\cdot A_{\\text{walls}} \\cdot (T_{\\text{mug}} - T_{\\text{room}})\\]"
                    </div>
                    <p>
                        "The mug's outer walls lose heat to the surrounding air. This is why "
                        "you can feel the warmth when holding a hot mug."
                    </p>
                </div>

                <div class="heat-flow">
                    <h3>"4. Evaporation from Coffee Surface"</h3>
                    <div class="formula">
                        "\\[Q_{\\text{evap}} = f_{\\text{open}} \\cdot \\eta \\cdot h_{\\text{mass}} \\cdot A_{\\text{top}} \\cdot L_v \\cdot (\\rho_{v,\\text{surface}} - \\rho_{v,\\text{air}})\\]"
                    </div>
                    <p>
                        "Water evaporating from the open coffee surface carries away latent heat "
                        "(\\(L_v \\approx 2.26\\) MJ/kg). The vapor density at the surface follows the "
                        "saturation curve (Tetens formula), while the ambient vapor density depends "
                        "on relative humidity (assumed 50%). The mass transfer coefficient "
                        "\\(h_{\\text{mass}}\\) is derived from the Lewis analogy, with an enhancement "
                        "factor \\(\\eta \\approx 2.5\\) to account for the buoyant thermal plume that "
                        "carries vapor away faster than a flat-plate model predicts."
                    </p>
                    <p>
                        "Evaporation is the " <strong>"dominant heat loss mechanism"</strong>
                        " for an open cup, responsible for roughly 40\u{2013}50% of total cooling. "
                        "This is why blowing on coffee works\u{2014}it accelerates evaporation."
                    </p>
                </div>

                <div class="heat-flow">
                    <h3>"5. Thermal Radiation"</h3>
                    <div class="formula">
                        "\\[Q_{\\text{rad,coffee}} = f_{\\text{open}} \\cdot \\varepsilon_{\\text{water}} \\, \\sigma \\, A_{\\text{top}} \\, (T_{\\text{coffee}}^4 - T_{\\text{room}}^4)\\]"
                    </div>
                    <div class="formula">
                        "\\[Q_{\\text{rad,mug}} = \\varepsilon_{\\text{mug}} \\, \\sigma \\, A_{\\text{outer}} \\, (T_{\\text{mug}}^4 - T_{\\text{room}}^4)\\]"
                    </div>
                    <p>
                        "Both the coffee surface and mug outer walls radiate heat according to the "
                        "Stefan\u{2013}Boltzmann law. Water has high emissivity (\\(\\varepsilon = 0.96\\)); "
                        "ceramic and glass are similar (\\(\\approx 0.90\\)), while stainless steel is "
                        "much lower (\\(\\approx 0.10\\)). Radiation from the coffee surface is scaled by "
                        "the lid opening fraction, while the mug radiates from its full outer surface "
                        "regardless. Radiation contributes roughly 15\u{2013}20% of total heat loss at "
                        "typical coffee temperatures."
                    </p>
                </div>
            </section>

            <section class="physics-detail">
                <h2>"Numerical Solution: Forward Euler"</h2>
                <p>
                    "One way to solve these coupled ODEs is to step forward in time numerically. "
                    "At each small time step " "\\(dt\\)" ", we compute the heat flows and update "
                    "the temperatures directly:"
                </p>
                <div class="formula">
                    "\\[\\Delta T = \\frac{Q \\cdot dt}{m \\cdot c}\\]"
                </div>
                <p class="formula-explanation">
                    "Where " "\\(mc\\)" " is the thermal mass (mass times specific heat capacity). "
                    "This represents how much energy is needed to change the temperature by 1 degree."
                </p>

                <h3>"Update Rules"</h3>
                <p>"At each time step, both temperatures are updated simultaneously:"</p>
                <div class="formula">
                    "\\[T_{\\text{coffee}}^{\\text{new}} = T_{\\text{coffee}} - \\frac{(Q_{\\text{to mug}} + Q_{\\text{conv}} + Q_{\\text{evap}} + Q_{\\text{rad}}) \\cdot dt}{C_{\\text{coffee}}}\\]"
                </div>
                <div class="formula">
                    "\\[T_{\\text{mug}}^{\\text{new}} = T_{\\text{mug}} + \\frac{(Q_{\\text{from coffee}} - Q_{\\text{conv}} - Q_{\\text{rad}}) \\cdot dt}{C_{\\text{mug}}}\\]"
                </div>
                <p>
                    "This is the " <strong>"forward Euler method"</strong>
                    "\u{2014}the simplest numerical ODE integrator. It approximates the continuous "
                    "derivatives by finite differences: " "\\(dT/dt \\approx \\Delta T / dt\\)" ". "
                    "Our implementation uses " "\\(dt = 0.6\\)" " seconds, which is small enough "
                    "for the lumped model to remain stable and accurate."
                </p>

                <h3>"Why Use a Numerical Method?"</h3>
                <p>
                    "With evaporation (exponential vapor pressure) and radiation (\\(T^4\\) dependence), "
                    "the system is nonlinear\u{2014}there is no exact closed-form solution. "
                    "The forward Euler approach handles this naturally:"
                </p>
                <ul>
                    <li>"It works for " <em>"any"</em> " ODE system, including nonlinear ones where "
                        "no closed-form solution exists"</li>
                    <li>"Adding new heat loss mechanisms (evaporation, radiation, a lid) just means "
                        "adding terms to the temperature update\u{2014}no rederivation needed"</li>
                    <li>"The same idea\u{2014}discretize space and time, step forward\u{2014}is how the "
                        "PDE model works too, just with a spatial grid added"</li>
                </ul>

                <div class="note">
                    <strong>"Foreshadowing: "</strong>
                    "When we move to the heat equation (PDE model), the spatial temperature distribution "
                    "means we can't reduce the problem to two ODEs. The numerical approach (finite differences) "
                    "becomes essential, and the analytical solution becomes far more involved."
                </div>
            </section>

            <section class="physics-detail">
                <h2>"Analytical Solution: Eigenvalue Decomposition"</h2>
                <p>
                    "The numerical model is nonlinear (evaporation and radiation), but we can "
                    <em>"linearize"</em> " these terms to obtain a closed-form solution. "
                    "Radiation is linearized as \\(Q_{\\text{rad}} \\approx h_{\\text{rad}} \\cdot A \\cdot (T - T_{\\text{room}})\\) "
                    "where \\(h_{\\text{rad}} = 4\\varepsilon\\sigma T_m^3\\). Evaporation is similarly "
                    "linearized by computing an effective heat transfer coefficient at the current "
                    "temperature. These linearized terms are folded into the convective coefficients, "
                    "yielding a 2\u{00D7}2 linear system with an exact solution."
                </p>
                <p>
                    "To keep the linearization accurate as the temperature drops, the solver uses "
                    <strong>"piecewise linearization"</strong>
                    ": every 12 seconds, it re-computes the linearized coefficients at the current "
                    "temperature, re-solves the eigenvalue problem, and continues from there. "
                    "This keeps the analytical and numerical solvers within ~0.1\u{00B0}C of each other "
                    "across the full simulation."
                </p>

                <h3>"Continuous ODEs"</h3>
                <p>
                    "Define deviations from room temperature: "
                    "\\(u = T_{\\text{coffee}} - T_{\\text{room}}\\)" ", "
                    "\\(v = T_{\\text{mug}} - T_{\\text{room}}\\)" ". "
                    "Then introduce four rate constants:"
                </p>
                <ul>
                    <li>"\\(a = k_{\\text{eff}} \\,/\\, C_{\\text{coffee}}\\)" " — coffee\u{2194}mug coupling"</li>
                    <li>"\\(b = h_{\\text{eff,coffee}} \\,/\\, C_{\\text{coffee}}\\)" " — coffee\u{2192}environment (convection + linearized radiation + linearized evaporation)"</li>
                    <li>"\\(c = k_{\\text{eff}} \\,/\\, C_{\\text{mug}}\\)" " — mug\u{2190}coffee coupling"</li>
                    <li>"\\(d = h_{\\text{eff,mug}} \\,/\\, C_{\\text{mug}}\\)" " — mug\u{2192}environment (convection + linearized radiation)"</li>
                </ul>
                <p>"The coupled system becomes:"</p>
                <div class="formula">
                    "\\[\\frac{du}{dt} = -(a+b)\\,u + a\\,v\\]"
                </div>
                <div class="formula">
                    "\\[\\frac{dv}{dt} = c\\,u - (c+d)\\,v\\]"
                </div>

                <h3>"Matrix Form and Eigenvalues"</h3>
                <p>"Writing this as a matrix equation:"</p>
                <div class="formula">
                    "\\[\\frac{d}{dt}\\begin{pmatrix} u \\\\ v \\end{pmatrix} = \\underbrace{\\begin{pmatrix} -(a+b) & a \\\\ c & -(c+d) \\end{pmatrix}}_{\\mathbf{M}} \\begin{pmatrix} u \\\\ v \\end{pmatrix}\\]"
                </div>
                <p>"The eigenvalues of " "\\(\\mathbf{M}\\)" " are:"</p>
                <div class="formula">
                    "\\[\\lambda_{1,2} = \\frac{1}{2}\\!\\left[-(a+b+c+d) \\pm \\sqrt{(a+b-c-d)^2 + 4ac}\\right]\\]"
                </div>
                <p>
                    "Both eigenvalues are real and negative (the discriminant is always non-negative "
                    "and the trace is always negative), so both modes decay exponentially."
                </p>

                <h3>"General Solution (per segment)"</h3>
                <p>
                    "With eigenvectors " "\\(\\mathbf{v}_1, \\mathbf{v}_2\\)" " the solution within each "
                    "segment is:"
                </p>
                <div class="formula">
                    "\\[\\begin{pmatrix} u(t) \\\\ v(t) \\end{pmatrix} = c_1 \\,\\mathbf{v}_1\\, e^{\\lambda_1 t} + c_2 \\,\\mathbf{v}_2\\, e^{\\lambda_2 t}\\]"
                </div>
                <p>
                    "The constants " "\\(c_1, c_2\\)" " are determined by the temperatures at the "
                    "start of each segment. At each segment boundary, the linearized coefficients "
                    "\\(b\\)" " and " "\\(d\\)" " are recomputed at the current temperature, giving "
                    "new eigenvalues and eigenvectors for the next segment."
                </p>

                <h3>"Two Physical Time Scales"</h3>
                <p>
                    "Within each segment, the two eigenvalues correspond to two distinct time scales "
                    "\\(\\tau_i = 1/|\\lambda_i|\\)" " (these shift from segment to segment "
                    "as the linearization point tracks the cooling temperature):"
                </p>
                <ul>
                    <li>
                        <strong>"Fast mode"</strong>
                        " (larger " "\\(|\\lambda|\\)" ", shorter " "\\(\\tau\\)" "): "
                        "coffee and mug equilibrating with each other. This mode decays quickly as "
                        "the two temperatures converge."
                    </li>
                    <li>
                        <strong>"Slow mode"</strong>
                        " (smaller " "\\(|\\lambda|\\)" ", longer " "\\(\\tau\\)" "): "
                        "the coupled coffee-mug system cooling together toward room temperature. "
                        "This is the dominant mode at later times."
                    </li>
                </ul>

                <div class="note">
                    <strong>"Key insight: "</strong>
                    "Preheating the mug reduces the initial temperature difference between coffee and mug, "
                    "which shrinks the amplitude of the fast mode. More of the system's energy starts "
                    "in the slow mode, which decays over a much longer time scale — keeping your "
                    "coffee warmer for longer."
                </div>
            </section>

            <section class="physics-detail">
                <h2>"Material Properties"</h2>
                <p>"Different materials have different thermal properties that affect cooling:"</p>

                <table class="materials-table">
                    <thead>
                        <tr>
                            <th>"Material"</th>
                            <th>"Specific Heat"<br/><span class="unit">"J/(kg K)"</span></th>
                            <th>"Density"<br/><span class="unit">"kg/m"<sup>"3"</sup></span></th>
                            <th>"Conductivity"<br/><span class="unit">"W/(m K)"</span></th>
                            <th>"Emissivity"</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td>"Ceramic"</td>
                            <td>"880"</td>
                            <td>"2400"</td>
                            <td>"1.0"</td>
                            <td>"0.90"</td>
                        </tr>
                        <tr>
                            <td>"Glass"</td>
                            <td>"840"</td>
                            <td>"2500"</td>
                            <td>"1.0"</td>
                            <td>"0.92"</td>
                        </tr>
                        <tr>
                            <td>"Stainless Steel"</td>
                            <td>"500"</td>
                            <td>"8000"</td>
                            <td>"16.0"</td>
                            <td>"0.10"</td>
                        </tr>
                        <tr>
                            <td>"Plastic"</td>
                            <td>"1200"</td>
                            <td>"1200"</td>
                            <td>"0.2"</td>
                            <td>"0.90"</td>
                        </tr>
                    </tbody>
                </table>

                <p class="note">
                    "Water (coffee) has a specific heat of 4186 J/(kg K), much higher than mug materials."
                </p>
            </section>

            <section class="physics-detail">
                <h2>"Thermal Mass Calculation"</h2>
                <p>"The mug is modeled as a cylinder where height = 2 × radius:"</p>

                <div class="formula">
                    "\\[V_{\\text{inner}} = \\pi r^2 h = 2\\pi r^3\\]"
                </div>
                <div class="formula">
                    "\\[r = \\left(\\frac{V}{2\\pi}\\right)^{1/3}\\]"
                </div>

                <p>"The mug's material volume is the outer cylinder minus the inner cylinder:"</p>
                <div class="formula">
                    "\\[V_{\\text{material}} = V_{\\text{outer}} - V_{\\text{inner}} + V_{\\text{bottom}}\\]"
                </div>

                <p>"Then thermal mass is:"</p>
                <div class="formula">
                    "\\[C_{\\text{thermal}} = V_{\\text{material}} \\cdot \\rho \\cdot c\\]"
                </div>
            </section>

            <section class="physics-detail">
                <h2>"Why Preheating Works"</h2>
                <p>
                    "When you preheat your mug, you're reducing the initial temperature difference "
                    "between the coffee and the mug. This means:"
                </p>
                <ul>
                    <li>"Less heat flows from coffee to mug initially (smaller temperature gradient)"</li>
                    <li>"The coffee retains more of its original heat energy"</li>
                    <li>"The effect is most pronounced in the first 5-10 minutes"</li>
                </ul>
                <p>
                    "A preheated mug also has stored thermal energy that helps buffer against heat loss, "
                    "keeping the coffee-mug system warmer for longer."
                </p>
            </section>

            <section class="physics-detail">
                <h2>"Model Limitations"</h2>
                <p>"Despite including five heat transfer mechanisms, the model still makes simplifying assumptions:"</p>
                <ul>
                    <li>"Uniform temperature within coffee and mug (lumped capacitance)"</li>
                    <li>"Lid modeled as a simple opening fraction " "\\(f_{\\text{open}}\\)" " that scales "
                        "top-surface losses (convection, evaporation, and radiation from the coffee surface)"</li>
                    <li>"Constant room temperature and 50% relative humidity"</li>
                    <li>"Simplified cylindrical geometry (height = 2 \u{00D7} radius)"</li>
                    <li>"Evaporation uses an empirical enhancement factor (\u{00D7}2.5) over the Lewis analogy"</li>
                    <li>"The analytical solver uses piecewise linearization (re-solved every 12s), "
                        "which stays within ~0.1\u{00B0}C of the numerical solver"</li>
                </ul>
            </section>

            <section class="interactive">
                <h2>"Try It Yourself"</h2>
                <p>"Adjust the parameters to see how they affect the cooling curves described above."</p>
                <LumpedPanel/>
            </section>

            <nav class="back-nav bottom">
                <A href="/theory">"< Back to Theory"</A>
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
