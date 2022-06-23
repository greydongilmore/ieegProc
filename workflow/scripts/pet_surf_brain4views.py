from nilearn import surface
from brain4views import plot_surf4

lh_over = surface.vol_to_surf(snakemake.input.pet, snakemake.input.lh_pial, radius=1, kind="ball")
rh_over = surface.vol_to_surf(snakemake.input.pet, snakemake.input.rh_pial, radius=1, kind="ball")

plot_surf4(
    [snakemake.input.lh_pial, snakemake.input.rh_pial],
    overlays=[lh_over, rh_over],
    avg_method="mean",
    colorbar=True,
    output_file=snakemake.output.pet_png,
	dpi=350
)
