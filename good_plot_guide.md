# Good Plot Guide

When designing a plot, there are several aspects to consider to ensure that it clearly conveys the relevant information and is also aesthetically pleasing. This guide is intended to provide a non-exhaustive checklist of things to consider when making a plot.

While some of the following advice is general, some is specific to plotting fields through two-dimensional filled contour plots.

Overall, figures should contain the minimum amount of information necessary for the figure to make sense with little or no external explanation.

## General
Here are some things to consider when making a figure:
1. Which plots or subplots are actually necessary?
2. Do different quantities need to be compared? Would plotting the difference be appropriate?
3. When plotting multiple quantities, would it be better to use subplots or multiple separate plots? *This is generally a trade-off between the amount of space available on a given figure against the necessity that subplots need to seen side-by-side. When comparing different quantities in different subplots or figures, try to use the same contours or scale between the subplots or figures.*

## Labels
- [ ] Should the figure have an overall title?
- [ ] If using subplots, each subplot should have a title. *When plotting fields, it is often helpful to add the minimum and maximum values of the field. This can be done quickly with tomplot using the `tomplot_field_title` routine.*
- [ ] Is the font used for text aesthetically pleasing? *To use latex within strings, pre-pended with `r`. Including the `set_tomplot_style` command quickly facilitates latex formatting.*
- [ ] Is the size of the text the appropriate size for the use of the plot? *When using tomplot, the can be specified using the argument to `set_tomplot_style`.*
- [ ] Are all labels and titles spaced appropriately from the axes?

## x and y axes
- [ ] Are the x and y axes labelled?
- [ ] Which ticks are used for the x and y axes? *Unless there is a good reason, consider only having ticks at the maximum and minimum values of the axes. This helps to declutter the plot.*
- [ ] Are tick labels shown to an appropriate number of significant figures?

## Subplots
- [ ] If contours and colours are shared between different subplots, can a single colour bar be used? *The tomplot `add_colorbar_fig` routine allows colour bars to be added to the figure without resizing individual subplots. In contrast, `add_colorbar_ax` adds a colour bar to a specific subplot.*
- [ ] Can subplots share axes? *By avoiding the need to label some axes or ticks, the figure can be decluttered.*
- [ ] Are subplots appropriate spaced from one another? *This can be altered using the `subplots_adjust` routine.

## Contours
- [ ] When plotting contours, are the values of the contours chosen to be nice, round values? *When trying to quickly create contours (possibly for unseen data), nice contours can be generated using the `tomplot_contours` routine.*
- [ ] Does the field have tiny oscillations around some contour? *It may be better to remove the contour for that value, using the `remove_contour` argument to the `tomplot_cmap` routine.

## Colours
- [ ] Are the colours used in the plot suitable for those with colour blindness? *Avoid red and green if possible. A plethora of colour maps can be found [here](https://matplotlib.org/stable/users/explain/colors/colormaps.html).*
- [ ] Have you considered whether your need plot should reproduce well in black and white? *The `tomplot_cmap` routine can generate colour maps that are rescaled to adjust the set of colours used. By rescaling colours to be less extreme, they may be clearer if converted to monochrome.*
- [ ] Do out-of-scale values need to be shown? *The `extend_cmap` argument to `tomplot_cmap` allows extrema to be shown in specific colours.*