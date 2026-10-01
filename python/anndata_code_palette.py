
# Replace some colors in the deault palette
replacement_colors_dict = {
    "Macrophage": "#ffffff",
    "T Cell": "#eeeeee",
}

obs_col = 'Annotations'
categories = pbmc.obs[obs_col].cat.categories
new_color_list = list(pbmc.uns[f'{obs_col}_colors'])

for cat, color in replacement_colors_dict.items():
    i = categories.get_loc(cat)
    new_color_list[i] = color

pbmc.uns[f'{obs_col}_colors'] = new_color_list
