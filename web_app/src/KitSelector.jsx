import { FormControl, InputLabel, MenuItem, Select } from '@mui/material'
import React from 'react'

function KitSelector({ kits, kitFilter, setKitFilter, ...formControlProps }) {
    return (
        <FormControl {...formControlProps}>
            <InputLabel>Kit</InputLabel>
            <Select
                label="Kit"
                value={kitFilter}
                onChange={(e) => setKitFilter(e.target.value)}
            >
                <MenuItem value="">None</MenuItem>
                {kits.map((kit) => <MenuItem key={kit} value={kit}>{kit}</MenuItem>)}
            </Select>
        </FormControl>
    )
}

export default KitSelector