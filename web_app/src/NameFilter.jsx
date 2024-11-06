import React from 'react'
import { FormControl, TextField } from '@mui/material'

function NameFilter({ nameFilter, setNameFilter, ...formControlProps }) {
    const [inputValue, setInputValue] = React.useState(nameFilter);

    React.useEffect(() => {
        const timeoutId = setTimeout(() => {
            setNameFilter(inputValue);
        }, 500);

        return () => clearTimeout(timeoutId);
    }, [inputValue, setNameFilter]);

    return (
        <FormControl {...formControlProps}>
            <TextField
                id="name-filter"
                variant="outlined"
                label="Filter by Name"
                placeholder="Enter name..."
                value={inputValue}
                onChange={(event) => setInputValue(event.target.value)}
            />
        </FormControl>
    )
}

export default NameFilter