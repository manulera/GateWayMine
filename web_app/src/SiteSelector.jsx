import React from 'react'
import { Autocomplete, TextField, Chip, FormControl } from '@mui/material'

function SiteSelector({ sites, setSelectedSites, ...formControlProps }) {
    return (
        <FormControl {...formControlProps}>
            <Autocomplete
                id="site-selector"
                multiple
                options={sites}
                onChange={(event, newValue) => {
                    setSelectedSites(newValue);
                }}
                renderInput={(params) => (
                    <TextField
                        {...params}
                        variant="outlined"
                        label="Select Sites"
                        placeholder="Sites"
                    />
                )}
                renderTags={(value, getTagProps) =>
                    value.map((option, index) => (
                        <Chip
                            label={option}
                            {...getTagProps({ index })}
                        />
                    ))
                }
            />
        </FormControl>
    )
}

export default SiteSelector