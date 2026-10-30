use super::constants::{H_WATER_SIDE, WATER_DENSITY, WATER_SPECIFIC_HEAT};

/// Material types for coffee mugs
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MugMaterial {
    Ceramic,
    Glass,
    StainlessSteel,
    Plastic,
}

impl MugMaterial {
    /// Specific heat capacity in J/(kg·K)
    pub fn specific_heat(&self) -> f64 {
        match self {
            MugMaterial::Ceramic => 880.0,
            MugMaterial::Glass => 840.0,
            MugMaterial::StainlessSteel => 500.0,
            MugMaterial::Plastic => 1200.0,
        }
    }

    /// Density in kg/m³
    pub fn density(&self) -> f64 {
        match self {
            MugMaterial::Ceramic => 2400.0,
            MugMaterial::Glass => 2500.0,
            MugMaterial::StainlessSteel => 8000.0,
            MugMaterial::Plastic => 1200.0,
        }
    }

    /// Thermal conductivity in W/(m·K)
    pub fn thermal_conductivity(&self) -> f64 {
        match self {
            MugMaterial::Ceramic => 1.0,
            MugMaterial::Glass => 1.0,
            MugMaterial::StainlessSteel => 16.0,
            MugMaterial::Plastic => 0.2,
        }
    }

    /// Emissivity (dimensionless, 0–1)
    pub fn emissivity(&self) -> f64 {
        match self {
            MugMaterial::Ceramic => 0.90,
            MugMaterial::Glass => 0.92,
            MugMaterial::StainlessSteel => 0.10,
            MugMaterial::Plastic => 0.90,
        }
    }

    /// Display name
    pub fn name(&self) -> &'static str {
        match self {
            MugMaterial::Ceramic => "Ceramic",
            MugMaterial::Glass => "Glass",
            MugMaterial::StainlessSteel => "Stainless Steel",
            MugMaterial::Plastic => "Plastic",
        }
    }
}

/// Parameters defining a mug and its environment
#[derive(Debug, Clone, PartialEq)]
pub struct MugParameters {
    pub material: MugMaterial,
    pub volume_ml: f64,
    pub wall_thickness_mm: f64,
    pub coffee_temp_c: f64,
    pub room_temp_c: f64,
    pub preheated: bool,
    pub preheat_temp_c: f64,
    /// Convective heat transfer coefficient in W/(m²·K).
    /// Typical range: 5–25 for natural convection in air.
    pub h_conv: f64,
    /// Surface emissivity (0–1). Defaults to material value but can be overridden.
    pub emissivity: f64,
    /// Lid openness (0.0–1.0): 0 = sealed lid, 1 = open mug with no lid.
    /// The mug is always present; exposure sets the open share of the top,
    /// which carries evaporation, convection and radiation from the coffee.
    pub exposure: f64,
}

impl MugParameters {
    /// Calculate the mass of the mug walls in kg
    pub fn mug_mass_kg(&self) -> f64 {
        // Approximate mug as a cylinder
        // volume_ml gives us the inner volume
        let inner_volume_m3 = self.volume_ml * 1e-6;

        // Assume height = 2 * radius for a typical mug shape
        // V = π * r² * h = π * r² * 2r = 2π * r³
        // r = (V / 2π)^(1/3)
        let inner_radius_m = (inner_volume_m3 / (2.0 * std::f64::consts::PI)).powf(1.0 / 3.0);
        let height_m = 2.0 * inner_radius_m;

        let wall_m = self.wall_thickness_mm * 1e-3;
        let outer_radius_m = inner_radius_m + wall_m;

        // Volume of ceramic = outer cylinder - inner cylinder + bottom
        let outer_volume = std::f64::consts::PI * outer_radius_m.powi(2) * (height_m + wall_m);
        let inner_volume = std::f64::consts::PI * inner_radius_m.powi(2) * height_m;
        let material_volume_m3 = outer_volume - inner_volume;

        material_volume_m3 * self.material.density()
    }

    /// Calculate the thermal mass of the mug in J/K
    pub fn mug_thermal_mass(&self) -> f64 {
        self.mug_mass_kg() * self.material.specific_heat()
    }

    /// Calculate the thermal mass of the coffee in J/K
    pub fn coffee_thermal_mass(&self) -> f64 {
        let coffee_mass_kg = self.volume_ml * 1e-6 * WATER_DENSITY;
        coffee_mass_kg * WATER_SPECIFIC_HEAT
    }

    /// Top (open coffee) surface area in m²
    pub fn top_surface_area_m2(&self) -> f64 {
        let inner_volume_m3 = self.volume_ml * 1e-6;
        let inner_radius_m = (inner_volume_m3 / (2.0 * std::f64::consts::PI)).powf(1.0 / 3.0);
        std::f64::consts::PI * inner_radius_m.powi(2)
    }

    /// Open coffee surface in m²: the lid opening's share of the top.
    pub fn open_top_area_m2(&self) -> f64 {
        self.exposure * self.top_surface_area_m2()
    }

    /// Inner wall area wetted by the coffee (side + base) in m², for the
    /// h = 2r cylinder: 2πr·2r + πr² = 5πr².
    pub fn wetted_area_m2(&self) -> f64 {
        let inner_volume_m3 = self.volume_ml * 1e-6;
        let inner_radius_m = (inner_volume_m3 / (2.0 * std::f64::consts::PI)).powf(1.0 / 3.0);
        let height_m = 2.0 * inner_radius_m;
        let side = 2.0 * std::f64::consts::PI * inner_radius_m * height_m;
        let base = std::f64::consts::PI * inner_radius_m.powi(2);
        side + base
    }

    /// Lumped coffee-to-mug conductance in W/K: the water-side film and half
    /// the wall's conduction resistance in series. Half, because the lumped mug
    /// temperature is the wall's mean, so conduction runs to the mid-plane.
    pub fn coffee_mug_conductance(&self) -> f64 {
        let area = self.wetted_area_m2();
        let wall_m = self.wall_thickness_mm * 1e-3;
        let r_film = 1.0 / (H_WATER_SIDE * area);
        let r_wall = (wall_m / 2.0) / (self.material.thermal_conductivity() * area);
        1.0 / (r_film + r_wall)
    }

    /// Outer mug surface area (side + bottom) in m²
    pub fn outer_surface_area_m2(&self) -> f64 {
        let inner_volume_m3 = self.volume_ml * 1e-6;
        let inner_radius_m = (inner_volume_m3 / (2.0 * std::f64::consts::PI)).powf(1.0 / 3.0);
        let height_m = 2.0 * inner_radius_m;
        let wall_m = self.wall_thickness_mm * 1e-3;
        let outer_radius_m = inner_radius_m + wall_m;

        let side = 2.0 * std::f64::consts::PI * outer_radius_m * height_m;
        let bottom = std::f64::consts::PI * outer_radius_m.powi(2);
        side + bottom
    }

    /// Total surface area for heat loss (approximate) in m²
    pub fn surface_area_m2(&self) -> f64 {
        self.top_surface_area_m2() + self.outer_surface_area_m2()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mug_mass_reasonable() {
        let params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let mass = params.mug_mass_kg();
        // A typical ceramic mug weighs 300-500g
        assert!(
            mass > 0.2 && mass < 0.8,
            "Mug mass {} kg seems unreasonable",
            mass
        );
    }

    fn params(material: MugMaterial, wall_mm: f64) -> MugParameters {
        MugParameters {
            material,
            volume_ml: 350.0,
            wall_thickness_mm: wall_mm,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: material.emissivity(),
            exposure: 0.2,
        }
    }

    #[test]
    fn conductance_scales_with_conductivity() {
        let g = |m| params(m, 5.0).coffee_mug_conductance();
        assert!(g(MugMaterial::StainlessSteel) > g(MugMaterial::Ceramic));
        assert!(g(MugMaterial::Ceramic) > g(MugMaterial::Plastic));
    }

    #[test]
    fn conductance_decreases_with_thickness() {
        let thin = params(MugMaterial::Ceramic, 2.0).coffee_mug_conductance();
        let thick = params(MugMaterial::Ceramic, 8.0).coffee_mug_conductance();
        assert!(thin > thick);
    }

    #[test]
    fn conductance_bounded_by_water_side() {
        let p = params(MugMaterial::StainlessSteel, 5.0);
        assert!(p.coffee_mug_conductance() < H_WATER_SIDE * p.wetted_area_m2());
    }

    #[test]
    fn conductance_matches_hand_computation() {
        // r = (350e-6 / 2π)^(1/3); A = 5πr²; G = 1/(1/(500·A) + 0.0025/(1·A))
        let r = (350e-6_f64 / (2.0 * std::f64::consts::PI)).powf(1.0 / 3.0);
        let area = 5.0 * std::f64::consts::PI * r * r;
        let expected = area / (1.0 / 500.0 + 0.0025 / 1.0);
        let p = params(MugMaterial::Ceramic, 5.0);
        assert!((p.wetted_area_m2() - area).abs() / area < 1e-9);
        assert!((p.coffee_mug_conductance() - expected).abs() / expected < 1e-9);
    }
}
