import "../component-stylesheets/Spinner.css"

export default function LoadingSpinner() {
    return (
      <div className="spinner-container">
        <div className="lds-roller"><div></div><div></div><div></div><div></div><div></div><div></div><div></div><div></div></div>
        <div className="spinner-label">Loading genes, please wait</div>
      </div>
    );
  }
