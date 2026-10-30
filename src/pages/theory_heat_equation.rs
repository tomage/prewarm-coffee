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
                    <li><strong>"Region I: Coffee"</strong>" (water) \u{2014} high thermal mass, moderate conductivity"</li>
                    <li><strong>"Region II: Mug"</strong>" (ceramic/glass) \u{2014} moderate conductivity, provides insulation"</li>
                    <li><strong>"Region III: Air"</strong>" \u{2014} very low conductivity, acts as thermal barrier. Heat conduction through the air layer is modeled, and convective and radiative loss at the mug\u{2019}s outer surface is a surface loss term there (Newton\u{2019}s law of cooling)."</li>
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
                    <h3>"Interfaces (r = r\u{1d62} and r = r\u{2092}): Continuity and Flux Balance"</h3>
                    <p>"At each material interface the temperature is continuous:"</p>
                    <div class="formula">
                        "\\[T_{\\text{left}} = T_{\\text{right}}\\]"
                    </div>
                    <p>
                        "The heat flux arriving from the inside either crosses into the next region or "
                        "leaves through the interface\u{2019}s surface loss:"
                    </p>
                    <div class="formula">
                        "\\[-k_1 \\frac{\\partial T}{\\partial r}\\bigg|_{r_i} = -k_2 \\frac{\\partial T}{\\partial r}\\bigg|_{r_i} + q_{\\text{top}}, \\qquad -k_2 \\frac{\\partial T}{\\partial r}\\bigg|_{r_o} = -k_3 \\frac{\\partial T}{\\partial r}\\bigg|_{r_o} + h_{\\text{eff}} \\left(T - T_{\\text{room}}\\right)\\]"
                    </div>
                    <p>
                        "At " "\\(r_i\\)" ", " "\\(q_{\\text{top}}\\)" " is the loss through the lid "
                        "opening (convection, radiation and evaporation from the open coffee surface), "
                        "spread over the sphere\u{2019}s inner surface. At " "\\(r_o\\)" ", "
                        "\\(h_{\\text{eff}} = h_{\\text{conv}} + 4\\varepsilon\\sigma T_m^3\\)" " is "
                        "convection plus linearized radiation from the mug\u{2019}s outer surface. "
                        "With no lid opening and no surface loss, both conditions reduce to plain flux "
                        "continuity, " "\\(k_{\\text{left}}\\, \\partial T/\\partial r = k_{\\text{right}}\\, \\partial T/\\partial r\\)" "."
                    </p>
                </div>

                <div class="heat-flow">
                    <h3>"Outer Boundary (r = r"<sub>"max"</sub>"): Room Temperature"</h3>
                    <div class="formula">
                        "\\[T(r_{\\max}, t) = T_{\\text{room}}\\]"
                    </div>
                    <p>
                        "The simulation uses a finite domain. Both solvers hold the air at "
                        "\\(r_{\\max}\\)" " at room temperature (a Dirichlet condition); the losses to "
                        "the room happen at the mug\u{2019}s outer surface, above. "
                        "The default " "\\(r_{\\max} = 10 \\times r_{\\text{outer}}\\)" " "
                        "closely approximates the true infinite-domain behavior; a smaller "
                        "\\(r_{\\max}\\)" " holds room-temperature air closer to the mug and cools it faster."
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
                <h2>"Numerical Solution: Finite Volumes"</h2>
                <p>
                    "The most practical approach is to discretize both space and time. We divide "
                    "the radial domain " "\\([0, r_{\\max}]\\)" " into grid points, fine near the "
                    "interfaces and coarse in the bulk. Each grid point owns a thin spherical "
                    "shell around it (its " <strong>"control volume"</strong> "), and the solver tracks "
                    "the heat flowing across the faces between neighbouring control volumes. "
                    "This " <strong>"finite-volume"</strong> " form conserves energy exactly: heat that "
                    "leaves one control volume enters the next."
                </p>
                <p>
                    "Time stepping uses the " <strong>"Crank\u{2013}Nicolson"</strong> " method\u{2014}"
                    "an implicit scheme that averages the spatial operator at the current and next time level:"
                </p>
                <div class="formula">
                    "\\[\\left(I - \\tfrac{\\Delta t}{2} L\\right) \\mathbf{T}^{n+1} = \\left(I + \\tfrac{\\Delta t}{2} L\\right) \\mathbf{T}^n + \\Delta t\\, \\mathbf{s}\\]"
                </div>
                <p>
                    "where " "\\(L\\)" " is the discretized diffusion operator in conservative form "
                    "\\(\\frac{1}{\\rho c \\, r^2} \\frac{\\partial}{\\partial r}\\!\\left(r^2 k \\frac{\\partial T}{\\partial r}\\right)\\)" " "
                    "and " "\\(\\mathbf{s}\\)" " holds the surface-loss source terms. "
                    "At each time step, this tridiagonal linear system is solved efficiently using the Thomas algorithm."
                </p>

                <h3>"Stability and Start-up"</h3>
                <p>
                    "Unlike the explicit (forward Euler) method, which requires "
                    "\\(\\Delta t < \\Delta r^2/(2\\alpha)\\)" ", "
                    "Crank\u{2013}Nicolson is " <strong>"unconditionally stable"</strong>
                    "\u{2014}there is no CFL constraint on the time step. It does let very stiff modes "
                    "flip sign from one step to the next instead of decaying, which shows up as ringing "
                    "at the mug wall right after the start. The solver therefore takes its first 4 steps "
                    "with backward Euler (" "\\(\\left(I - \\Delta t\\, L\\right) \\mathbf{T}^{n+1} = \\mathbf{T}^n + \\Delta t\\, \\mathbf{s}\\)" "), "
                    "which damps those modes, and then switches to Crank\u{2013}Nicolson "
                    "(a " <strong>"Rannacher start-up"</strong> ")."
                </p>
                <p>
                    "The time step is " "\\(\\Delta t = \\min(t_{\\text{end}}/500,\\ 0.5\\,\\text{s})\\)" ", "
                    "and never longer than the interval between recorded snapshots. For any run longer "
                    "than about 4 minutes that is 0.5 s: 3,600 steps for 30 minutes, 14,400 for 2 hours."
                </p>

                <h3>"Interface Handling"</h3>
                <p>
                    "At material boundaries, " "\\(\\alpha\\)" " and " "\\(k\\)" " change abruptly. "
                    "A grid point sits exactly on each interface, so every grid segment lies in one "
                    "material and the flux across it uses that material\u{2019}s " "\\(k\\)" ". "
                    "The interface point\u{2019}s control volume is part one material and part the other, "
                    "and its heat capacity is the sum of the two parts. The heat arriving through one "
                    "face and leaving through the other is balanced in that shared volume, which enforces "
                    "the flux balance at " "\\(r_i\\)" " and " "\\(r_o\\)" " without approximating "
                    "\\(k\\)" " across the jump."
                </p>

                <h3>"Boundary Conditions"</h3>
                <p>
                    "The grid point at " "\\(r_{\\max}\\)" " is held at " "\\(T_{\\text{room}}\\)" ". "
                    "At the shell\u{2019}s outer surface, convection and radiation (linearized at the "
                    "current surface temperature, re-done every step) act on the interface point. At the "
                    "coffee\u{2019}s surface " "\\(r_i\\)" ", the open-top loss is linearized every step "
                    "as the tangent " "\\(h\\,(T - T_{\\text{ref}})\\)" " at the current temperature. "
                    "Evaporation keeps going at room temperature, so " "\\(T_{\\text{ref}}\\)" " is a little "
                    "below room and the coffee can settle slightly below room temperature."
                </p>

                <div class="note">
                    <strong>"Comparison with lumped model: "</strong>
                    "The lumped model has 2 unknowns (coffee temp, mug temp) and takes microseconds. "
                    "The PDE solver has hundreds of grid points and takes milliseconds\u{2014}still fast enough "
                    "for interactive use, but capturing the spatial structure that the lumped model misses."
                </div>

                <div class="note">
                    <strong>"Why start with numerical? "</strong>
                    "Finite volumes handle material discontinuities naturally, "
                    "and adding nonlinear effects like radiation, evaporation, and convective "
                    "surface losses is straightforward. "
                    "The analytical eigenfunction expansion (below) provides exact results for "
                    "the linear problem but cannot handle these nonlinearities directly; it "
                    "linearizes them afresh for each of ten time segments (see below)\u{2014}"
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
                    "For a single homogeneous sphere whose surface is held at "
                    "\\(T_{\\text{room}}\\)" ", the eigenvalues are simply "
                    "\\(\\lambda_n = n\\pi / R\\)" " and the full solution is:"
                </p>
                <div class="formula">
                    "\\[T(r,t) = T_{\\text{room}} + \\sum_{n=1}^{\\infty} A_n \\frac{\\sin(n\\pi r/R)}{n\\pi r/R} \\, e^{-\\alpha (n\\pi/R)^2 t}\\]"
                </div>
                <p>
                    "The coefficients " "\\(A_n\\)" " are determined by the initial excess over room, "
                    "\\(f(r) = T(r,0) - T_{\\text{room}}\\)" ", using orthogonality of the eigenfunctions "
                    "(with weight " "\\(r^2\\)" "):"
                </p>
                <div class="formula">
                    "\\[A_n = \\frac{2 n\\pi}{R^2} \\int_0^R r \\, f(r) \\, \\sin\\!\\left(\\frac{n\\pi r}{R}\\right) dr\\]"
                </div>
                <p>
                    "For a uniform initial temperature " "\\(T(r,0) = T_0\\)" ", this gives "
                    "\\(A_n = 2 (T_0 - T_{\\text{room}}) (-1)^{n+1}\\)" "."
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
                    "\\(g_n = -\\sin(\\beta_3 r_{\\max})\\)" " make " "\\(T_3 = T_{\\text{room}}\\)" " at "
                    "\\(r_{\\max}\\)" ". (On an infinite domain, the spectrum becomes continuous and "
                    "\\(T_3\\)" " would involve an integral transform rather than a discrete sum.)"
                </p>
                <p>
                    "The eigenvalues are no longer simple multiples of " "\\(\\pi\\)" ". Instead, "
                    "matching temperature continuity and the flux balance at each interface gives four "
                    "conditions relating " "\\(A_n, B_n, C_n, D_n\\)" ". For non-trivial solutions, the 4\u{00D7}4 "
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
                <p>
                    "An eigenfunction expansion needs every coefficient to be constant, but the surface "
                    "losses change as the coffee cools. The analytical solver therefore cuts the run into "
                    <strong>"ten equal time segments"</strong> " and holds the losses constant within each. "
                    "At " "\\(r_i\\)" ", the open-top loss becomes "
                    "\\(q_{\\text{top}} = h_{\\text{inner}} (T - T_{\\text{room}})\\)" ", with "
                    "\\(h_{\\text{inner}}\\)" " the open share of the sphere times the loss at the coffee\u{2019}s "
                    "average temperature divided by its excess over room (taken at 5 K above room if the "
                    "coffee is closer to room than that). At " "\\(r_o\\)" ", " "\\(h_{\\text{eff}}\\)"
                    " takes the radiation term at the mug\u{2019}s outer-surface temperature."
                </p>
                <p>
                    "Each segment uses the " <strong>"midpoint rule"</strong> ": it is solved once with the "
                    "losses taken at the segment\u{2019}s starting temperatures, to predict where it ends, and "
                    "solved again with the losses taken halfway between start and predicted end. Each "
                    "segment is exact for its constant coefficients. The next one starts from its final "
                    "temperature profile, re-expanded in the new segment\u{2019}s eigenfunctions by projecting "
                    "the profile onto them (Simpson\u{2019}s rule in each region). The numerical solver "
                    "re-linearizes every step."
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
                        " If T"<sub>"l"</sub>" > T"<sub>"s"</sub>", heat flows from coffee \u{2192} mug \u{2192} air. "
                        "If T"<sub>"l"</sub>" < T"<sub>"s"</sub>", the mug initially heats the coffee."
                    </li>
                    <li>
                        <strong>"Water\u{2019}s thermal mass:"</strong>
                        " Water has high specific heat (4186 J/kg\u{00B7}K) \u{2014} it stores lots of heat but transfers it slowly."
                    </li>
                    <li>
                        <strong>"Ceramic as insulator:"</strong>
                        " The ceramic mug has low thermal conductivity, slowing heat transfer to the air."
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
                <h2>"Assumptions"</h2>
                <ul>
                    <li>"The mug is a sphere: the coffee is a sphere of the chosen volume, inside a "
                        "spherical shell of the chosen wall thickness, surrounded by air."</li>
                    <li>"Still air with a constant convection coefficient " "\\(h_{\\text{conv}}\\)" " at the "
                        "mug\u{2019}s outer surface and the open coffee surface."</li>
                    <li>"Water\u{2019}s effective conductivity is 6.0 W/(m\u{00B7}K), about 10\u{00D7} "
                        "molecular, standing in for the convective mixing inside the coffee."</li>
                    <li>"Radiation is linearized: " "\\(h_{\\text{rad}} = 4\\varepsilon\\sigma T_m^3\\)" ". The "
                        "numerical solver re-linearizes every step at the current surface temperature; the "
                        "analytical solver re-linearizes once per segment, at the outer surface\u{2019}s "
                        "midpoint temperature for that segment."</li>
                    <li>"Evaporation follows the Lewis analogy \u{00D7} 2.5 for the rising plume, at 50% "
                        "relative humidity. It continues at room temperature, so in the numerical solver "
                        "coffee can end a few degrees below room. The 2.5 factor scales evaporation only, not "
                        "convection, which puts that limit slightly below the wet-bulb temperature."</li>
                    <li>"A lid is a perfect insulator with no heat capacity: no heat leaves through the covered "
                        "part of the coffee\u{2019}s surface."</li>
                    <li>"At 100% exposure (no lid), the open area is a typical mug\u{2019}s top, mapped to "
                        "19.1% of the sphere\u{2019}s inner surface; exposure scales it linearly. The open-top "
                        "losses are spread evenly over the inner surface, because the model is radially symmetric."</li>
                    <li>"The analytical solver linearizes the open-top losses per segment, at the coffee\u{2019}s "
                        "average temperature for that segment, but at least 5 K above room. It can therefore "
                        "only cool toward room temperature, never below it."</li>
                    <li>"The analytical solver takes the open-top loss at the coffee\u{2019}s average "
                        "temperature; the numerical one at the open surface, which is cooler. So the analytical "
                        "coffee runs slightly cold: its temperature differs from the numerical one\u{2019}s by "
                        "up to about 1.8 \u{00B0}C with no lid (180 ml, " "\\(h_{\\text{conv}}\\)" " = 20) and "
                        "by about 0.2 \u{00B0}C or less with the lid at most 20% open. For a ceramic mug at the "
                        "default 20% opening (180 or 500 ml; 2, 5 or 8 mm walls) the two agree on what "
                        "preheating buys to within 0.1 minute of drinkable time and 0.2 \u{00B0}C, but the "
                        "analytical coffee drops below drinkable 0.2\u{2013}0.4 minutes early (larger mugs "
                        "more), in both the cold and the preheated mug. The gap widens with the lid half or "
                        "fully open, where the crossings land about 0.4 minutes early on average and up to "
                        "1 minute, and for stainless steel at low convection, where the analytical solver "
                        "credits preheating with up to 0.5 minute more than the numerical one."</li>
                    <li>"The air domain is finite: air at " "\\(r_{\\max}\\)" " (default 10\u{00D7} the mug\u{2019}s "
                        "outer radius) is held at room temperature."</li>
                    <li>"Constant material properties, independent of temperature; constant room temperature."</li>
                </ul>
            </section>

            <section class="physics-detail">
                <h2>"Model Limitations"</h2>
                <p>"Compared to the lumped model, the PDE approach removes the uniform-temperature "
                   "assumption but has limits of its own:"</p>
                <ul>
                    <li>"Radial symmetry only (no vertical convection currents, no distinct top surface)"</li>
                    <li>"Natural convection in the liquid is approximated via an effective conductivity, not resolved with fluid dynamics"</li>
                    <li>"Pure conduction within the air region (convective loss is modeled at the mug\u{2019}s outer surface, but not internal air circulation)"</li>
                </ul>
            </section>

            <section class="interactive">
                <h2>"Try It Yourself"</h2>
                <p>"Adjust the parameters to see how heat flows through the coffee, ceramic mug, and surrounding air."</p>
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
