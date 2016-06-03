First, copy the TCU input bits (layer2 output bits) from the trigger page's
DSM algorithm pdf files (trigger page -> trigger system logic on left side)

For each DSM, make a dat file with columns `[output bit] [trigger name]`

Then make a `tcu.dat` which has columns `[dsm] [output bit] [trigger name]`

Finally, remove "Unused" and undesired bits from `tcu.dat`

Now create a file `tcuchan.dat` which has columns `[DSM] [tcu input channel]`,
read from `L1_Crate_Cable_Map.pdf` found on the same trigger system logic page

The `TCUbits` class ONLY reads data from `tcu.dat` and `tcuchan.dat`
