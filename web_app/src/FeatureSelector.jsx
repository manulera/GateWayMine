import React from 'react'
import { Autocomplete, TextField, Chip, FormControl } from '@mui/material'

function FeatureSelector({ features, setSelectedFeatures, ...formControlProps }) {
    return (
        <FormControl {...formControlProps}>
            <Autocomplete
                id="feature-selector"
                multiple
                options={features.filter(feature => !(/^att[BPLR]\d$/.test(feature)))}
                onChange={(event, newValue) => {
                    setSelectedFeatures(newValue);
                }}
                renderInput={(params) => (
                    <TextField
                        {...params}
                        variant="outlined"
                        label="Select Features"
                        placeholder="Features"
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

export default React.memo(FeatureSelector);