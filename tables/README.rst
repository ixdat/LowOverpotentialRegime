This is the "relational database" at the center of the repository. Each folder should
be thought of as a table and each .json file should be thought of as a row in that table.
The reason for doing it this way is to keep the database human-readable without any
programming required. The .json files are text files which can be opened with notepad.

This may be superseded with an ixdat sqlite database in time (ixdat v0.3.x)