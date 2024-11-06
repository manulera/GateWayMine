import React from 'react'
import SiteSelector from './SiteSelector'
import FeatureSelector from './FeatureSelector'
import { Box } from '@mui/material'
import NameFilter from './NameFilter'
import SourceSelector from './SourceSelector.jsx'
import KitSelector from './KitSelector'

function FilterBar({ sites, features, setSelectedSites, setSelectedFeatures, nameFilter, setNameFilter, selectedSource, setSelectedSource, kits, kitFilter, setKitFilter }) {
    return (
        <Box sx={{ mx: 'auto', display: 'flex', flexWrap: 'wrap', gap: 2, justifyContent: 'center' }}>
            <NameFilter nameFilter={nameFilter} setNameFilter={setNameFilter} sx={{ width: 200 }} />
            <SiteSelector sites={sites} setSelectedSites={setSelectedSites} sx={{ width: 200 }} />
            <FeatureSelector features={features} setSelectedFeatures={setSelectedFeatures} sx={{ width: 200 }} />
            <SourceSelector selectedSource={selectedSource} setSelectedSource={setSelectedSource} sx={{ width: 200 }} />
            <KitSelector kits={kits} kitFilter={kitFilter} setKitFilter={setKitFilter} sx={{ width: 200 }} />
        </Box>
    )
}

export default FilterBar