use crate::components::HeatEquationPanel;
use leptos::*;
use leptos_router::A;

#[component]
pub fn TheoryHeatEquationPage() -> impl IntoView {
    view! {
        <main class="container details-page">
            <header>
                <h1>"Theory: Heat Equation (PDE)"</h1>
                <p class="tagline">"Consider a "<a href="https://en.wikipedia.org/wiki/Spherical_cow" target="_blank" rel="noopener">"spherical"</a>" coffee cup in space\u{2026}"</p>
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
                <h2>"The Problem"</h2>
                <p>
                    "The lumped capacitance model assumes uniform temperature within each body. "
                    "In reality, the center of your coffee is hotter than the edges, and temperature "
                    "varies continuously through the mug wall. To capture these spatial gradients, "
                    "we need the heat equation\u{2014}a partial differential equation (PDE)."
                </p>
                <p>
                    "We model the system as a sphere with three concentric regions:"
                </p>
                <ul>
                    <li><strong>"Region I: Liquid"</strong>" (water/coffee) \u{2014} high thermal mass, moderate conductivity"</li>
                    <li><strong>"Region II: Shell"</strong>" (ceramic/glass) \u{2014} moderate conductivity, provides insulation"</li>
                    <li><strong>"Region III: Air"</strong>" \u{2014} very low conductivity, acts as thermal barrier. Heat conduction through the air layer is modeled, and convective loss at the outer surface is captured via a Robin boundary condition (Newton\u{2019}s law of cooling)."</li>
                </ul>
            </section>

            <section class="physics-detail">
                <h2>"The Governing Equation"</h2>
                <p>"In spherical coordinates with radial symmetry, the heat equation is:"</p>
                <div class="formula">
                    "\\[\\frac{\\partial T}{\\partial t} = \\alpha \\left( \\frac{\\partial^2 T}{\\partial r^2} + \\frac{2}{r} \\frac{\\partial T}{\\partial r} \\right)\\]"
                </div>
                <p class="formula-explanation">
                    "Where " "\\(\\alpha = k / (\\rho c)\\)" " is the thermal diffusivity\u{2014}a measure of "
                    "how quickly heat spreads through a material. Each region has its own " "\\(\\alpha\\)" "."
                </p>
                <p>
                    "Unlike the lumped model (which produces ODEs in time only), this PDE involves "
                    "derivatives in both space (" "\\(r\\)" ") and time (" "\\(t\\)" "). "
                    "The temperature " "\\(T(r, t)\\)" " is now a function of position and time."
                </p>
            </section>

            <section class="physics-detail">
                <h2>"Boundary Conditions"</h2>
                <p>"The PDE needs conditions at each boundary and interface:"</p>

                <div class="heat-flow">
                    <h3>"Center (r = 0): Symmetry"</h3>
                    <div class="formula">
                        "\\[\\frac{\\partial T}{\\partial r}\\bigg|_{r=0} = 0\\]"
                    </div>
                    <p>"No heat flows through the center of the sphere."</p>
                </div>

                <div class="heat-flow">
                    <h3>"Interfaces (r = r\u{1d62} and r = r\u{2092}): Continuity"</h3>
                    <p>"At each material interface, both temperature and heat flux must be continuous:"</p>
                    <div class="formula">
                        "\\[T_{\\text{left}} = T_{\\text{right}}, \\quad k_{\\text{left}} \\frac{\\partial T}{\\partial r} = k_{\\text{right}} \\frac{\\partial T}{\\partial r}\\]"
                    </div>
                    <p>"The flux condition accounts for different thermal conductivities on each side."</p>
                </div>

                <div class="heat-flow">
                    <h3>"Outer Boundary (r = r"<sub>"max"</sub>"): Robin Condition"</h3>
                    <div class="formula">
                        "\\[-k_3 \\frac{\\partial T}{\\partial r}\\bigg|_{r_{\\max}} = h_{\\text{eff}} \\left(T - T_{\\text{room}}\\right)\\]"
                    </div>
                    <p>
                        "The simulation uses a finite domain. At the outer boundary, a Robin condition "
                        "models convective and radiative loss to the environment. With "
                        "\\(r_{\\max} \\geq 10 \\times r_{\\text{outer}}\\)" ", "
                        "this closely approximates the true infinite-domain behavior."
                    </p>
                </div>
            </section>

            <section class="physics-detail">
                <h2>"Material Properties"</h2>
                <p>"The three regions have very different thermal properties:"</p>

                <table class="materials-table">
                    <thead>
                        <tr>
                            <th>"Material"</th>
                            <th>"k"<br/><span class="unit">"W/(m\u{00B7}K)"</span></th>
                            <th>"\u{03C1}"<br/><span class="unit">"kg/m\u{00B3}"</span></th>
                            <th>"c"<br/><span class="unit">"J/(kg\u{00B7}K)"</span></th>
                            <th>"\u{03B1}"<br/><span class="unit">"m\u{00B2}/s"</span></th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td>"Water"</td>
                            <td>"6.0*"</td>
                            <td>"1000"</td>
                            <td>"4186"</td>
                            <td>"1.43\u{00D7}10\u{207B}\u{2076}"</td>
                        </tr>
                        <tr>
                            <td>"Ceramic"</td>
                            <td>"1.0"</td>
                            <td>"2400"</td>
                            <td>"880"</td>
                            <td>"4.73\u{00D7}10\u{207B}\u{2077}"</td>
                        </tr>
                        <tr>
                            <td>"Air"</td>
                            <td>"0.026"</td>
                            <td>"1.18"</td>
                            <td>"1005"</td>
                            <td>"2.19\u{00D7}10\u{207B}\u{2075}"</td>
                        </tr>
                    </tbody>
                </table>

                <p class="note">
                    "*" <strong>"Effective conductivity:"</strong>
                    " The molecular conductivity of water is only 0.6 W/(m\u{00B7}K), "
                    "but hot coffee has vigorous buoyancy-driven convection (Rayleigh number ~ 10\u{2079}) "
                    "that keeps the bulk nearly well-mixed. The simulation uses an effective conductivity "
                    "of ~10\u{00D7} molecular to approximate this mixing without solving the full "
                    "Navier\u{2013}Stokes equations. This also raises the effective thermal diffusivity by 10\u{00D7}."
                </p>
                <p class="note">
                    "Air has high diffusivity but extremely low conductivity, "
                    "making it an excellent insulator."
                </p>
            </section>

            <section class="physics-detail">
                <h2>"Numerical Solution: Finite Differences"</h2>
                <p>
                    "The most practical approach is to discretize both space and time. We divide "
                    "the radial domain into grid points on a finite domain " "\\([0, r_{\\max}]\\)" " and "
                    "step forward in time using the " <strong>"Crank\u{2013}Nicolson"</strong> " method\u{2014}"
                    "an implicit scheme that averages the spatial operator at the current and next time level:"
                </p>
                <div class="formula">
                    "\\[\\left(I - \\tfrac{\\Delta t}{2} L\\right) \\mathbf{T}^{n+1} = \\left(I + \\tfrac{\\Delta t}{2} L\\right) \\mathbf{T}^n + \\Delta t\\, \\mathbf{s}\\]"
                </div>
                <p>
                    "where " "\\(L\\)" " is the discretized diffusion operator in conservative form "
                    "\\(\\frac{1}{\\rho c \\, r^2} \\frac{\\partial}{\\partial r}\\!\\left(r^2 k \\frac{\\partial T}{\\partial r}\\right)\\)" " "
                    "and " "\\(\\mathbf{s}\\)" " includes convection and radiation source terms. "
                    "At each time step, this tridiagonal linear system is solved efficiently using the Thomas algorithm."
                </p>

                <h3>"Unconditional Stability"</h3>
                <p>
                    "Unlike the explicit (forward Euler) method, which requires "
                    "\\(\\Delta t < \\Delta r^2/(2\\alpha)\\)" ", "
                    "Crank\u{2013}Nicolson is " <strong>"unconditionally stable"</strong>
                    "\u{2014}there is no CFL constraint on the time step. "
                    "This allows the solver to use ~500 time steps regardless of grid resolution, "
                    "keeping computation fast even with fine spatial grids."
                </p>

                <h3>"Interface Handling"</h3>
                <p>
                    "At material boundaries, " "\\(\\alpha\\)" " and " "\\(k\\)" " change abruptly. "
                    "The solver uses the " <strong>"harmonic mean"</strong> " of conductivities at half-points "
                    "\\(k_{i+1/2} = 2 k_i k_{i+1} / (k_i + k_{i+1})\\)" " to properly enforce "
                    "flux continuity " "\\(k_1 \\partial T/\\partial r = k_2 \\partial T/\\partial r\\)" " at interfaces."
                </p>

                <h3>"Boundary Conditions"</h3>
                <p>
                    "The simulation uses a finite domain with a Dirichlet condition "
                    "\\(T(r_{\\max}) = T_{\\text{room}}\\)" " at the outer boundary. "
                    "Convection and linearized radiation are applied as Robin conditions "
                    "at the shell outer surface. With " "\\(r_{\\max} \\geq 10 \\times r_{\\text{outer}}\\)" ", "
                    "the finite domain has negligible effect on liquid/shell temperatures."
                </p>

                <div class="note">
                    <strong>"Comparison with lumped model: "</strong>
                    "The lumped model has 2 unknowns (coffee temp, mug temp) and takes microseconds. "
                    "The PDE solver has hundreds of grid points and takes milliseconds\u{2014}still fast enough "
                    "for interactive use, but capturing the spatial structure that the lumped model misses."
                </div>

                <div class="note">
                    <strong>"Why start with numerical? "</strong>
                    "Finite differences handle material discontinuities naturally, "
                    "and adding nonlinear effects like radiation, evaporation, and convective "
                    "surface losses is straightforward. "
                    "The analytical eigenfunction expansion (below) provides exact results for "
                    "the linear problem but cannot handle these nonlinearities directly\u{2014}"
                    "both approaches are available in this project."
                </div>
            </section>

            <section class="physics-detail">
                <h2>"Analytical Solution: Eigenfunction Expansion"</h2>
                <p>
                    "The heat equation also has an exact solution, but deriving it is considerably "
                    "more involved than the 2\u{00D7}2 eigenvalue problem of the lumped model."
                </p>

                <h3>"The Substitution Trick"</h3>
                <p>
                    "The key insight is that the substitution " "\\(u(r,t) = r \\cdot T(r,t)\\)" " transforms the "
                    "spherical heat equation into the standard 1D form. Starting with " "\\(T = u/r\\)" ":"
                </p>
                <div class="formula">
                    "\\[\\frac{\\partial T}{\\partial t} = \\frac{1}{r}\\frac{\\partial u}{\\partial t}, \\quad \\frac{\\partial T}{\\partial r} = \\frac{1}{r}\\frac{\\partial u}{\\partial r} - \\frac{u}{r^2}, \\quad \\frac{\\partial^2 T}{\\partial r^2} = \\frac{1}{r}\\frac{\\partial^2 u}{\\partial r^2} - \\frac{2}{r^2}\\frac{\\partial u}{\\partial r} + \\frac{2u}{r^3}\\]"
                </div>
                <p>
                    "Substituting into the heat equation and simplifying, the cross terms cancel, leaving:"
                </p>
                <div class="formula">
                    "\\[\\frac{\\partial u}{\\partial t} = \\alpha \\frac{\\partial^2 u}{\\partial r^2}\\]"
                </div>

                <h3>"Separation of Variables"</h3>
                <p>
                    "Assuming " "\\(u(r,t) = R(r) \\cdot \\Theta(t)\\)" " separates the PDE into:"
                </p>
                <ul>
                    <li><strong>"Time part:"</strong>" " "\\(\\Theta(t) = e^{-\\alpha \\lambda^2 t}\\)" " \u{2014} exponential decay"</li>
                    <li><strong>"Spatial part:"</strong>" " "\\(R(r) = A \\sin(\\lambda r) + B \\cos(\\lambda r)\\)" " \u{2014} oscillatory"</li>
                </ul>
                <p>
                    "Regularity at " "\\(r = 0\\)" " forces " "\\(B = 0\\)" " (since "
                    "\\(\\cos(\\lambda r)/r \\to \\infty\\)" "), giving eigenfunctions of the form "
                    "\\(\\sin(\\lambda r) / r\\)" "\u{2014}the zeroth-order spherical Bessel function."
                </p>

                <h3>"Single Region: Clean Eigenvalues"</h3>
                <p>
                    "For a single homogeneous sphere, the eigenvalues are simply "
                    "\\(\\lambda_n = n\\pi / R\\)" " and the full solution is:"
                </p>
                <div class="formula">
                    "\\[T(r,t) = T_{\\text{room}} + \\sum_{n=1}^{\\infty} A_n \\frac{\\sin(n\\pi r/R)}{n\\pi r/R} \\, e^{-\\alpha (n\\pi/R)^2 t}\\]"
                </div>
                <p>
                    "The coefficients " "\\(A_n\\)" " are determined by the initial condition "
                    "\\(T(r,0) = f(r)\\)" " using orthogonality of the eigenfunctions:"
                </p>
                <div class="formula">
                    "\\[A_n = \\frac{2}{R} \\int_0^R r \\, f(r) \\, \\sin\\!\\left(\\frac{n\\pi r}{R}\\right) dr\\]"
                </div>
                <p>
                    "For a uniform initial temperature " "\\(T(r,0) = T_0\\)" ", this gives "
                    "\\(A_n = 2 T_0 (-1)^{n+1} / (n\\pi)\\)" "."
                </p>

                <h3>"Three Regions: Transcendental Eigenvalues"</h3>
                <p>
                    "With three regions, each having different " "\\(\\alpha\\)" " and " "\\(k\\)" ", "
                    "the solution in each region has its own form."
                </p>
                <p>
                    <strong>"Notation change: "</strong>
                    "In the single-region case, " "\\(\\lambda\\)" " is a spatial frequency (units 1/m) and the "
                    "time decay is " "\\(e^{-\\alpha \\lambda^2 t}\\)" ". With three regions this is awkward\u{2014}"
                    "each region has a different " "\\(\\alpha\\)" ", yet all regions must share the same "
                    "decay rate for a given eigenmode. The standard approach is to redefine "
                    "\\(\\lambda_n\\)" " so that " "\\(\\lambda_n^2\\)" " already has units of 1/s, "
                    "with each region\u{2019}s " "\\(\\alpha\\)" " absorbed into its spatial wave number: "
                    "region I uses " "\\(\\lambda_n / \\sqrt{\\alpha_1}\\)" " as its spatial frequency, "
                    "region II uses " "\\(\\lambda_n / \\sqrt{\\alpha_2}\\)" ", and so on. "
                    "The time part then becomes simply " "\\(e^{-\\lambda_n^2 t}\\)" " everywhere."
                </p>
                <div class="formula">
                    "\\[T_1(r,t) = T_{\\text{room}} + \\sum_n A_n \\frac{\\sin(\\lambda_n r / \\sqrt{\\alpha_1})}{r} \\, e^{-\\lambda_n^2 t}\\]"
                </div>
                <div class="formula">
                    "\\[T_2(r,t) = T_{\\text{room}} + \\sum_n \\frac{B_n \\sin(\\cdots) + C_n \\cos(\\cdots)}{r} \\, e^{-\\lambda_n^2 t}\\]"
                </div>
                <div class="formula">
                    "\\[T_3(r,t) = T_{\\text{room}} + \\sum_n D_n \\frac{g_d \\sin(\\beta_3 r) + g_n \\cos(\\beta_3 r)}{r} \\, e^{-\\lambda_n^2 t}\\]"
                </div>
                <p class="formula-explanation">
                    "where " "\\(\\beta_3 = \\lambda_n / \\sqrt{\\alpha_3}\\)" ". "
                    "On a finite domain, region III uses standing waves (sin/cos) rather than "
                    "exponential decay. The coefficients "
                    "\\(g_d = \\cos(\\beta_3 r_{\\max})\\)" " and "
                    "\\(g_n = -\\sin(\\beta_3 r_{\\max})\\)" " enforce the boundary condition at "
                    "\\(r_{\\max}\\)" ". (On an infinite domain, the spectrum becomes continuous and "
                    "\\(T_3\\)" " would involve an integral transform rather than a discrete sum.)"
                </p>
                <p>
                    "The eigenvalues are no longer simple multiples of " "\\(\\pi\\)" ". Instead, "
                    "matching temperature and flux continuity at each interface gives four conditions "
                    "relating " "\\(A_n, B_n, C_n, D_n\\)" ". For non-trivial solutions, the 4\u{00D7}4 "
                    "determinant must vanish:"
                </p>
                <div class="formula">
                    "\\[\\det\\!\\left(\\mathbf{M}(\\lambda)\\right) = 0\\]"
                </div>
                <p>
                    "This is a " <strong>"transcendental equation"</strong> "\u{2014}it cannot be solved "
                    "in closed form. The eigenvalues must be found numerically (bisection, Newton's method), "
                    "and the series must be truncated at a finite number of terms."
                </p>

                <div class="note">
                    <strong>"Key contrast with the lumped model: "</strong>
                    "The lumped model's 2\u{00D7}2 eigenvalue problem has an explicit formula for "
                    "\\(\\lambda_{1,2}\\)" ". The PDE's eigenvalue problem produces infinitely many "
                    "eigenvalues, each requiring a numerical root-find. The analytical solution exists "
                    "in principle, but in practice still requires numerical computation."
                </div>
            </section>

            <section class="physics-detail">
                <h2>"Physical Interpretation"</h2>

                <h3>"Decay Modes"</h3>
                <p>
                    "Each eigenvalue " "\\(\\lambda_n\\)" " corresponds to a spatial mode that decays "
                    "with time constant " "\\(\\tau_n = 1/\\lambda_n^2\\)" " (recall that in the three-region "
                    "convention, " "\\(\\lambda_n^2\\)" " already has units of 1/s):"
                </p>
                <ul>
                    <li>
                        <strong>"Fundamental mode"</strong>" (" "\\(\\lambda_1\\)" ", slowest decay): "
                        "represents the overall cooling of the system. Dominates at long times."
                    </li>
                    <li>
                        <strong>"Higher modes"</strong>" (" "\\(\\lambda_n\\)" " for " "\\(n > 1\\)" "): "
                        "represent spatial gradients (center vs. edge). These decay quickly, "
                        "which is why temperature profiles flatten out early."
                    </li>
                </ul>

                <h3>"Characteristic Time Scale"</h3>
                <p>
                    "The fundamental mode gives the classic scaling:"
                </p>
                <div class="formula">
                    "\\[\\tau \\sim \\frac{R^2}{\\pi^2 \\alpha}\\]"
                </div>
                <p>
                    "Larger objects (bigger " "\\(R\\)" ") and materials with lower diffusivity "
                    "(smaller " "\\(\\alpha\\)" ") take longer to cool. This explains why a large "
                    "mug of coffee stays warm much longer than a thin espresso cup."
                </p>
            </section>

            <section class="physics-detail">
                <h2>"Key Observations"</h2>
                <ul>
                    <li>
                        <strong>"Heat flow direction:"</strong>
                        " If T"<sub>"l"</sub>" > T"<sub>"s"</sub>", heat flows from liquid \u{2192} shell \u{2192} air. "
                        "If T"<sub>"l"</sub>" < T"<sub>"s"</sub>", the shell initially heats the liquid."
                    </li>
                    <li>
                        <strong>"Water\u{2019}s thermal mass:"</strong>
                        " Water has high specific heat (4186 J/kg\u{00B7}K) \u{2014} it stores lots of heat but transfers it slowly."
                    </li>
                    <li>
                        <strong>"Ceramic as insulator:"</strong>
                        " The ceramic shell has low thermal conductivity, slowing heat transfer to the air."
                    </li>
                    <li>
                        <strong>"Air as insulator:"</strong>
                        " Air\u{2019}s low conductivity creates a thermal barrier \u{2014} this is why real systems use convection."
                    </li>
                    <li>
                        <strong>"Preheating effect:"</strong>
                        " Try setting T"<sub>"s"</sub>" close to T"<sub>"l"</sub>" \u{2014} the liquid stays hot longer "
                        "(same physics as preheating your coffee mug!)."
                    </li>
                </ul>
            </section>

            <section class="physics-detail">
                <h2>"Model Limitations"</h2>
                <p>"Compared to the lumped model, the PDE approach removes the uniform-temperature "
                   "assumption but introduces others:"</p>
                <ul>
                    <li>"Spherical geometry (real mugs are cylindrical)"</li>
                    <li>"Radial symmetry only (no vertical convection currents)"</li>
                    <li>"Natural convection in the liquid is approximated via an effective conductivity (10\u{00D7} molecular), not resolved with fluid dynamics"</li>
                    <li>"Pure conduction within the air region (convective loss is modeled at the outer boundary, but not internal air circulation)"</li>
                    <li>"Evaporative cooling is included in the numerical solver (controlled by the exposure slider) but handled approximately in the analytical solver via a separate exposed-sphere model"</li>
                    <li>"Finite simulation domain (not truly infinite)\u{2014}mitigated by using " "\\(r_{\\max} \\geq 10 \\times r_{\\text{outer}}\\)"</li>
                </ul>
            </section>

            <section class="interactive">
                <h2>"Try It Yourself"</h2>
                <p>"Adjust the parameters to see how heat flows through the liquid, ceramic shell, and surrounding air."</p>
                <HeatEquationPanel/>
            </section>

            <section class="physics-detail">
                <h2>"Further Reading"</h2>
                <ul>
                    <li>"Carslaw & Jaeger, " <em>"Conduction of Heat in Solids"</em> " \u{2014} the classic reference"</li>
                    <li>"Incropera & DeWitt, " <em>"Fundamentals of Heat and Mass Transfer"</em></li>
                    <li>"Haberman, " <em>"Applied Partial Differential Equations"</em> " \u{2014} good for mathematical details"</li>
                </ul>
            </section>

            <nav class="back-nav bottom">
                <A href="/theory">"< Back to Theory"</A>
            </nav>

            <footer>
                <p>
                    <a href="https://github.com/tomage/prewarm-coffee" target="_blank" rel="noopener">"Source code"</a>
                    " · AGPL-3.0 · "
                    "2025\u{2013}2026"
                </p>
            </footer>
        </main>
    }
}
