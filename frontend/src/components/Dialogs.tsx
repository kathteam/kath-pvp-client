import {
  Dialog,
  DialogActions,
  DialogTitle,
  DialogContent,
  DialogContentText,
  TextField,
  Button,
} from '@mui/material';

export function RenameDialog({ open, onClose, onConfirm }: { open: boolean; onClose: () => void; onConfirm: () => void }) {
  return (
    <Dialog open={open} onClose={onClose}>
      <DialogTitle>Rename File</DialogTitle>
      <DialogContent>
        <DialogContentText>Enter a new name for the file:</DialogContentText>
        <TextField autoFocus margin="dense" id="new-name" label="New Name" type="text" fullWidth variant="standard" />
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Cancel</Button>
        <Button onClick={onConfirm}>Rename</Button>
      </DialogActions>
    </Dialog>
  );
}

export function DeleteDialog({ open, onClose, onConfirm }: { open: boolean; onClose: () => void; onConfirm: () => void }) {
  return (
    <Dialog open={open} onClose={onClose}>
      <DialogTitle>Delete File</DialogTitle>
      <DialogContent>
        <DialogContentText>Are you sure you want to delete the file?</DialogContentText>
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Cancel</Button>
        <Button onClick={onConfirm}>Delete</Button>
      </DialogActions>
    </Dialog>
  );
}

export function CreateDatabaseDialog({
  open,
  onClose,
  dbFilename,
  setDbFilename,
  onConfirm,
}: {
  open: boolean;
  onClose: () => void;
  dbFilename: string;
  setDbFilename: (value: string) => void;
  onConfirm: () => void;
}) {
  return (
    <Dialog open={open} onClose={onClose}>
      <DialogTitle>Create Database</DialogTitle>
      <DialogContent>
        <DialogContentText>Enter a name for the database file:</DialogContentText>
        <TextField
          autoFocus
          margin="dense"
          id="db-filename"
          label="Database Filename"
          type="text"
          fullWidth
          variant="standard"
          value={dbFilename}
          onChange={(e) => setDbFilename(e.target.value)}
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Cancel</Button>
        <Button onClick={onConfirm}>Create</Button>
      </DialogActions>
    </Dialog>
  );
}