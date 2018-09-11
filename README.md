# CEE484_HW4_MSE_wall_design

Description:
    The script was written to primarily for design of MSE walls.
    The script calculates minimum length of geogrid required to resist loads from the pressure from 
    (1) earth and (2) surcharge with groundwater (GW) level is low to not effect the MSE wall.
Assumption and Considerations:
    (1) groundwater (GW) level is very low - does not effect the MSE wall stability.
    (2) global stability is adequate.
    (3) the surcharge load is uniform surcharge; no point load or line load
    (4) geogrid is assumed to be extandable; hence, follows failure wedge defined by Rankine theory
        for internal stability check. 
    (5) all the geogrid length and vertical spacing is kept constant for simplicity.
    (6) pull-out resistance factor of geogrid is assumed to be 0.7*tan(effective friction angle of MSE fill)
    (7) inclination of wall and backfill is vertical and horizontal, respectively.
    (8) minimum width-to-height ratio (B/H) for MSE wall is 0.7
