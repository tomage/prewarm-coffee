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
#[derive(Debug, Clone)]
pub struct MugParameters {
    pub material: MugMaterial,
    pub volume_ml: f64,
    pub wall_thickness_mm: f64,
    pub coffee_temp_c: f64,
    pub room_temp_c: f64,
    pub preheated: bool,
    pub preheat_temp_c: f64,
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
        // Water: specific heat = 4186 J/(kg·K), density ≈ 1 kg/L
        let coffee_mass_kg = self.volume_ml * 1e-3; // 1 ml water ≈ 1 g
        coffee_mass_kg * 4186.0
    }

    /// Surface area for heat loss (approximate) in m²
    pub fn surface_area_m2(&self) -> f64 {
        let inner_volume_m3 = self.volume_ml * 1e-6;
        let inner_radius_m = (inner_volume_m3 / (2.0 * std::f64::consts::PI)).powf(1.0 / 3.0);
        let height_m = 2.0 * inner_radius_m;
        let wall_m = self.wall_thickness_mm * 1e-3;
        let outer_radius_m = inner_radius_m + wall_m;

        // Side surface + top (open) + bottom
        let side = 2.0 * std::f64::consts::PI * outer_radius_m * height_m;
        let top = std::f64::consts::PI * inner_radius_m.powi(2); // coffee surface
        let bottom = std::f64::consts::PI * outer_radius_m.powi(2);

        side + top + bottom
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
        };

        let mass = params.mug_mass_kg();
        // A typical ceramic mug weighs 300-500g
        assert!(
            mass > 0.2 && mass < 0.8,
            "Mug mass {} kg seems unreasonable",
            mass
        );
    }
}
