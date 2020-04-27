# Mengsong-wood-decompo-FEM
These are the files that contain data used to analyze and to reproduce Dossa et al manuscript if Forest Ecology and Management. These are from an over 3 years wood decomposition experiment across a tropical forest disturbance landscape in Mengsong, SW, China. 
These are in total 8 files.
1	Mengsong_field_data_FEM.csv which is the main file and consists of data collected from the two native species logs used in the experiment. 
1.1	Tag_n.x unique identifier of wood log,
1.2	For_type forest type with OC = open canopy (regenerating forest), CC = closed canopy (mature forest) and OL = open land,
1.3	PLOT experimental plot unit,
1.4	S-PLOT sub unit embedded within the experimental plot,
1.5	tree_index the identifier of individual tree were logs come from. Numbers for Litsea cubeba and letters for Castanopsis mekongensis,
1.6	Species_full full latin name of wood species used,
1.7	Species, Species.x short name of wood species used,
1.8	Pos_soil relative position to the soil where wood core was collected from. Up = half core not in direct contact with soil and down = half wood core in direct contact with soil,
1.9	Bark_thickness_T1 through T4 bark thickness (mm) measured with electronic caliper two measurements at each end of the log,
1.10	Installment_date date when wood log was incubated on the forest floor,
1.11	Thickness_bark_fresh1.mm_later and Thickness_bark_fresh2.mm_later thickness (mm) of fresh bark at the end of the experiment. One measurement taken at each end of log,
1.12	Bark_Dry_weight initial disk dried weight (g),
1.13	Disk_Fresh_weight_g_ini initial disk fresh weight (g),
1.14	Disk_Green_volume_ini initial disk fresh volume (g) measured by water displacement method (Williamson & Wiemann, 2010),
1.15	Disk_Dry_weight_ini initial disk dried weight (g),
1.16	Disk_Dry_volume_ini initial disk dried volume (g) measured by water displacement method (Williamson & Wiemann, 2010),
1.17	Water_cont_ini.y intial log water content calculated as (fresh weight-dry weight)/dry weight (Jones et al., 2019),
1.18	WSG_ini.y wood specific gravity calculated as (oven dry weight/oven dry volume)* water density see (Williamson & Wiemann, 2010),
1.19	Bark_sample_Fresh.weight.1_ini initial disk bark sample fresh weight (g),
1.20	Bark_sample_Green_volume.ini initial disk bark sample fresh volume(g) measured by water displacement method (Williamson & Wiemann, 2010),
1.21	Bark_sample_Dry_weight.ini initial disk bark sample dried weight (g),
1.22	Bark_sample_Dry_volume.ini initial disk bark sample dried volume (g) measured by water displacement method (Williamson & Wiemann, 2010),
1.23	db, dm and dt diameter of wood log a bottom (end 1), middle, and top (end 2) respectively (cm),
1.24	L wood log length (cm),
1.25	Coll_No.x collection or harvest sequence. 1 for first harvest at 3 months, and 6 for sixth harvest at 36 months. In between harvest happened at 6 month interval, no harvest was done at 30 months,
1.26	Coll_date.x collection or harvest date,
1.27	Number_days number of days spanned since wood log incubation (days),
1.28	Fresh_mass_later harvested wood core fresh weight (g),
1.29	Fresh_volume_later harvested wood core fresh volume (g) measured by water displacement method (Williamson & Wiemann, 2010),
1.30	Oven_dry_mass.105.for.60h_later harvested wood core dried mass (g),
1.31	Oven_dry_volume_later harvested wood core dried volume (g) measured by water displacement method (Williamson & Wiemann, 2010),
1.32	Water_cont_later.y harvested wood core water content calculated as (fresh weight-dry weight)/dry weight (Jones et al., 2019),
1.33	WSG_later.y harvested wood core wood specific gravity calculated as (oven dry weight/oven dry volume)* water density see (Williamson & Wiemann, 2010), 
1.34	Per_WSG_loss percentage of wood specific gravity loss calculated as (wood specific gravity initial- wood specific gravity later)*100/wood specific gravity initial,
1.35	Log_volume wood wood log volume (cm^3) calculated with Netwon formula (Harmon & Sexton, 1996),
1.36	Wood_density wood density (g cm-3) calculated as dried weight/fresh volume (Mori et al., 2013),
1.37	Log_mass_ini log initial dried mass (g) calculated knowing wood density and volume,
1.38	Log_Bark_Dry_weight_g log remaining bark dried weight (g) at the end of experiment,
1.39	Log_No_Bark_Dry_weight_g log remaining dried weight (g) with bark removed,
1.40	Log_mass_final remaining log dried weight (g) at the end of the experiment, calculated as remaining log dried weight + log remaining bark dried weight at the end of experiment,
1.41	Bark_sample_Fresh.weight_later sample of remaining bark fresh weight (g) at the end of the experiment,
1.42	Bark_sample_Green.volume_later sample of remaining bark fresh volume (g) at the end of the experiment measured by water displacement method (Williamson & Wiemann, 2010),
1.43	Bark_sample_Dry.weight_later sample of remaining bark dried weight (g) at the end of experiment,
1.44	Bark_sample_Dry.volume_later sample of remaining bark dried volume (g) at the end of experiment measured by water displacement method (Williamson & Wiemann, 2010),
1.45	ML mass loss (g),
1.46	ML_percent percentage mass loss which is calculated as mass loss * 100/initial mass,
1.47	m_initial log initial dried weight (g),
1.48	m_harvest log dried weight (g) at the end of the experiment,
1.49	t number of days (days) spanned on forest floor till harvest,
1.50	Termi.assum log termite status presence (1), absence (0) with assumption that termite stay after first encounter or record in the log during fieldwork,
1.51	PLOT_SPLOT combination of PLOT unit with S_PLOT sub plot unit. This was used later to combine this data file with the environmental data file. 
2	Mengsong_environmental_factors.csv contains the wood samples chemistry analysis content in percent (%).
3	Preddat.csv represents the new dataframe to be used to predict wood specific gravity loss based on selected best statistical model, 
4	Termite_infestation_per_plot.csv compiles the number of logs infested by termites (n_presence), the number of logs non-infested by termites (n_absence) per wood species (Species.x), per harvest (Coll_No.x.x) and per plot (PLOT.x), 
5	Environmental factors_Mengsong.csv which comprises all the soil/topography and vegetation of plots and subplots. Soil/topography data were reduced to the first 2 PCA axes with original analyses in (Paudel et al., 2015). Vegetation data consists of the first 3 NMDS axes with original analyses in (Paudel et al., 2015),
6	Temperature_landscape_government station.csv compares the microclimate temperature (recorded with Hobo microclimate station) with the local climate recorded by local government station,
7	Temp_Landscape_only.csv consists of microclimate temperature (recorded by Hobo microclimate station) installed in each of the disturbance gradient category,
8	Soil moisture and RH_Landscape.csv compiles both daily average soil moisture and ambient relative humidity (recorded by Hobo microclimate station) installed in each of the disturbance gradient category,
9	PAR_median_plot.csv consists of the photosynthetically active radiation (recorded by Hobo microclimate station) for daily median of readings taken 1 hour either side of the solar noon.  


References cited
Harmon, M. E., & Sexton, J. (1996). Guidelines for Measurements of Woody Detritus in Forest Ecosystems. Publication No. 20, U.S.LTER Network Office: University of Washington, Seattle, WA, USA.
Jones, J. M., Heath, K. D., Ferrer, A., Brown, S. P., Canam, T., & James, W. (2019). Wood decomposition in aquatic and terrestrial ecosystems in the tropics: contrasting biotic and abiotic processes. FEMS Microbiology Ecology, 95, fiy223.
Mori, S., Itoh,  a., Nanami, S., Tan, S., Chong, L., & Yamakura, T. (2013). Effect of wood density and water permeability on wood decomposition rates of 32 Bornean rainforest trees. Journal of Plant Ecology, 1–8. doi: 10.1093/jpe/rtt041
Paudel, E., Dossa, G. G. O., Blécourt, M. de, Beckschäfer, P., Xu, J., Harrison, R. D., … Harrison, R. D. (2015). Quantifying the factors affecting leaf litter decomposition across a tropical forest disturbance gradient. Ecosphere, 6(December), 267.
Williamson, G. B., & Wiemann, M. C. (2010). Measuring wood specific gravity...Correctly. American Journal of Botany, 97(3), 519–524. doi: 10.3732/ajb.0900243

