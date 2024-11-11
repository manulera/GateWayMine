import React from 'react'
import { FormControl, InputLabel, Select, MenuItem } from '@mui/material'

function SourceSelector({ selectedSource, setSelectedSource, ...formControlProps }) {
    return (
        <FormControl {...formControlProps}>
            <InputLabel id="source-select-label">Source</InputLabel>
            <Select
                labelId="source-select-label"
                id="source-select"
                value={selectedSource}
                label="Source"
                onChange={(e) => setSelectedSource(e.target.value)}
            >
                <MenuItem value="all">All Sources</MenuItem>
                <MenuItem value="addgene">Addgene</MenuItem>
                <MenuItem value="snapgene">Snapgene</MenuItem>
            </Select>
        </FormControl>
    )
}

export default React.memo(SourceSelector);